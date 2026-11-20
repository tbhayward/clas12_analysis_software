"""Check the supplied theory-column mapping and asymmetric bands."""
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plot_kroll_comparisons as plot


class PredictionChecks(unittest.TestCase):
    def setUp(self):
        self.predictions = plot.load_new_theory(plot.DEFAULT_PREDICTION_CSV)

    def test_supplied_ratios_are_not_rescaled_or_sign_flipped(self):
        raw = pd.read_csv(plot.DEFAULT_PREDICTION_CSV, comment='#').sort_values(['xB','minus_tprime']).reset_index(drop=True)
        for observable, column in plot.NEW_THEORY_COLUMNS.items():
            for suffix in ('', '_lo', '_hi'):
                np.testing.assert_array_equal(self.predictions[observable+suffix], raw[column+suffix])
        np.testing.assert_array_equal(self.predictions.groupby('xB').Q2.first(), [2.2,2.6,3.1,4.3])
        self.assertLess(self.predictions.loc[0, 'lu1'], 0)

    def test_band_uses_error_lengths_and_only_supplied_abscissae(self):
        pred = self.predictions[self.predictions.xB == .2]
        fig, ax = plt.subplots()
        try:
            plot.draw_new_theory_curve(ax, pred, 'll1', 'black')
            np.testing.assert_array_equal(ax.lines[0].get_xdata(), pred.minus_tprime)
            np.testing.assert_array_equal(ax.lines[0].get_ydata(), pred.ll1)
            vertices = ax.collections[0].get_paths()[0].vertices
            for _, row in pred.iterrows():
                ys = vertices[np.isclose(vertices[:,0],row.minus_tprime),1]
                self.assertAlmostEqual(ys.min(), row.ll1-row.ll1_lo)
                self.assertAlmostEqual(ys.max(), row.ll1+row.ll1_hi)
        finally:
            plt.close(fig)

    def test_malformed_inputs_are_rejected(self):
        raw = pd.read_csv(plot.DEFAULT_PREDICTION_CSV, comment='#')
        cases = {'duplicate': pd.concat([raw,raw.iloc[:1]],ignore_index=True),
                 'missing_column':raw.drop(columns=['A_LL_hi']),
                 'missing_unpolarized':raw.drop(columns=['A_UU_cosphi'])}
        negative=raw.copy();negative.loc[0,'A_LL_lo']=-.01;cases['negative_error']=negative
        nonfinite=raw.copy();nonfinite.loc[0,'A_LL']=np.nan;cases['nonfinite']=nonfinite
        mixed=raw.copy();mixed.loc[0,'Q2']=9.;cases['multiple_q2']=mixed
        with tempfile.TemporaryDirectory() as directory:
            for name,frame in cases.items():
                with self.subTest(name=name):
                    path=Path(directory)/(name+'.csv');frame.to_csv(path,index=False)
                    with self.assertRaises(ValueError):plot.load_new_theory(path)

    def test_axes_contain_the_supplied_bands(self):
        for observable in plot.NEW_THEORY_OBSERVABLES:
            low,high=plot.new_theory_y_limits(self.predictions,observable)
            self.assertLessEqual(low,(self.predictions[observable]-self.predictions[observable+'_lo']).min())
            self.assertGreaterEqual(high,(self.predictions[observable]+self.predictions[observable+'_hi']).max())

    def test_measurement_error_versions(self):
        from unittest.mock import Mock
        frame = pd.DataFrame({'mean_minus_tprime_gev2':[.2], 'll1':[.1],
                              'll1_stat':[.03], 'll1_point_to_point_systematic':[.04]})
        combined = Mock()
        plot.draw_measurements(combined,frame,'ll1','black','o',combined_errors=True)
        np.testing.assert_allclose(combined.errorbar.call_args.kwargs['yerr'],[.05])
        combined.add_patch.assert_not_called()
        original = Mock()
        plot.draw_measurements(original,frame,'ll1','black','o')
        np.testing.assert_allclose(original.errorbar.call_args.kwargs['yerr'],[.03])
        self.assertEqual(original.add_patch.call_count,1)

    def test_no_combined_theory_all_xb_canvas(self):
        with self.assertRaises(ValueError):
            plot.plot_new_theory_canvas(None,None,Path('.'),kroll_predictions=plot.kroll_frame())


if __name__=='__main__':
    unittest.main()
