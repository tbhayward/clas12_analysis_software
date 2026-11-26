"""Full-BM charged-pion comparisons for the two target-spin observables."""
from pathlib import Path
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.pyplot as plt


def draw(theory,data,output,dpi,alpha,obs,combined_errors=False):
 import plot_comparisons as p
 from figure_layout import OBS_NAMES
 fig,axes=plt.subplots(2,2,figsize=(11.7,8.7))
 for xi,ax in enumerate(axes.flat):
  p.axis(ax,obs,(-.4,.4) if obs=='ul1' else (-1,1))
  lo,hi=p.XB[xi];ax.set_title(rf'${lo:.2f}<x_B<{hi:.2f}$')
  for sp,c,m,ls in [('pip','#0072B2','o','-'),('pim','#D55E00','s','--')]:
   b=theory[sp,obs][xi]['arrays'];ax.plot(b['minus_tprime'],b['total'],color=c,ls=ls,lw=2.2)
   p.points(ax,p.state_data(data,sp,obs,xi),'combined_errors' if combined_errors else 'stat_syst' if sp=='pip' else 'stat',c,m)
 handles=[Line2D([],[],color=c,marker=m,mfc='white',ls=ls,label=p.CHANNEL[sp]+' data / full BM') for sp,c,m,ls in [('pip','#0072B2','o','-'),('pim','#D55E00','s','--')]]
 if not combined_errors:handles.append(Patch(facecolor='#0072B2',alpha=.15,label=r'$\pi^+$ point-to-point syst.'))
 fig.legend(handles=handles,loc='lower center',ncol=2,frameon=False)
 fig.suptitle(p.LABELS[obs]+' · Charged-pion comparison',fontsize=18)
 fig.subplots_adjust(left=.095,right=.98,bottom=.16,top=.9,hspace=.32,wspace=.23);p.watermark(fig,alpha)
 folder=Path(output)/'figures/charge_comparison';folder.mkdir(parents=True,exist_ok=True)
 path=folder/(f'pi_plus_pi_minus_{OBS_NAMES[obs]}'+('_combined_errors' if combined_errors else '')+'_2x2.png')
 fig.savefig(path,dpi=dpi);plt.close(fig)
 return dict(path=str(path.relative_to(output)),title='π+ / π− '+OBS_NAMES[obs],category='charge comparison')
