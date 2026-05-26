#!/usr/bin/env python3

import argparse
import math
import re
import sys
from collections import OrderedDict


def fatal(message):
    print(f"[FATAL] {message}", file=sys.stderr)
    sys.exit(1)


def parse_triples(rhs_text, variable_name):
    """
    Parse a Mathematica-style list of triples:

      {{x1, y1, dy1}, {x2, y2, dy2}, ...}

    into:

      [(x1, y1, dy1), (x2, y2, dy2), ...]
    """

    triple_pattern = re.compile(
        r"\{\s*"
        r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*,\s*"
        r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*,\s*"
        r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*"
        r"\}"
    )

    triples = []
    for match in triple_pattern.finditer(rhs_text):
        x_val = float(match.group(1))
        y_val = float(match.group(2))
        dy_val = float(match.group(3))

        triples.append((x_val, y_val, dy_val))

    #endfor

    if len(triples) == 0:
        fatal(f"Could not parse any triples for variable '{variable_name}'.")

    #endif

    return triples


def parse_file(path):
    """
    Parse all assignments of the form:

      variableName = {{...}, {...}, ...};

    Returns an OrderedDict:

      variable_name -> list of triples
    """

    try:
        with open(path, "r") as input_file:
            text = input_file.read()

    except OSError as exc:
        fatal(f"Could not read input file '{path}': {exc}")

    #endtry

    assignment_pattern = re.compile(
        r"([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(\{\{.*?\}\})\s*;",
        re.DOTALL
    )

    results = OrderedDict()

    for match in assignment_pattern.finditer(text):
        variable_name = match.group(1)
        rhs_text = match.group(2)

        if variable_name in results:
            fatal(f"Duplicate assignment for variable '{variable_name}' in file '{path}'.")

        #endif

        results[variable_name] = parse_triples(rhs_text, variable_name)

    #endfor

    if len(results) == 0:
        fatal(f"No valid assignments found in file '{path}'.")

    #endif

    return results


def validate_same_structure(file_data, input_paths):
    """
    Enforce that all files contain exactly the same variable names,
    in the same order, and with the same number of triples.
    """

    reference_names = list(file_data[0].keys())

    for file_index in range(1, len(file_data)):
        current_names = list(file_data[file_index].keys())

        if current_names != reference_names:
            missing = sorted(set(reference_names) - set(current_names))
            extra = sorted(set(current_names) - set(reference_names))

            message_parts = [
                f"File '{input_paths[file_index]}' does not match the variable list in '{input_paths[0]}'."
            ]

            if missing:
                message_parts.append(f"Missing variables: {missing}")

            #endif

            if extra:
                message_parts.append(f"Extra variables: {extra}")

            #endif

            fatal(" ".join(message_parts))

        #endif

    #endfor

    for variable_name in reference_names:
        reference_length = len(file_data[0][variable_name])

        for file_index in range(1, len(file_data)):
            current_length = len(file_data[file_index][variable_name])

            if current_length != reference_length:
                fatal(
                    f"Variable '{variable_name}' has {current_length} triples in "
                    f"'{input_paths[file_index]}', but {reference_length} triples in "
                    f"'{input_paths[0]}'."
                )

            #endif

        #endfor

    #endfor


def weighted_average_triplet(triples, variable_name, point_index):
    """
    Combine one corresponding triple from each input file.

    Each input triple is:

      (mean_tprime, asymmetry_value, uncertainty)

    The weighted average uses weights:

      w_i = 1 / sigma_i^2

    for both mean_tprime and asymmetry_value.

    The propagated uncertainty is:

      sigma_avg = 1 / sqrt(sum_i w_i)

    Special case:
      If every uncertainty is exactly zero, this script requires all values
      to be identical. It then returns that identical value with zero uncertainty.
    """

    zero_uncertainty_count = sum(1 for _, _, dy_val in triples if dy_val == 0.0)

    if zero_uncertainty_count == len(triples):
        first_x = triples[0][0]
        first_y = triples[0][1]

        for x_val, y_val, dy_val in triples:
            if x_val != first_x or y_val != first_y:
                fatal(
                    f"All uncertainties are zero for variable '{variable_name}', "
                    f"point index {point_index}, but the central values are not identical. "
                    f"Cannot compute a finite weighted average."
                )

            #endif

        #endfor

        return first_x, first_y, 0.0

    #endif

    if zero_uncertainty_count > 0:
        fatal(
            f"Variable '{variable_name}', point index {point_index}, has a mixture of "
            f"zero and nonzero uncertainties. Refusing to silently choose a convention."
        )

    #endif

    sum_weights = 0.0
    weighted_x_sum = 0.0
    weighted_y_sum = 0.0

    for x_val, y_val, dy_val in triples:
        if dy_val < 0.0:
            fatal(
                f"Variable '{variable_name}', point index {point_index}, has a negative "
                f"uncertainty: {dy_val}."
            )

        #endif

        weight = 1.0 / (dy_val * dy_val)

        sum_weights += weight
        weighted_x_sum += weight * x_val
        weighted_y_sum += weight * y_val

    #endfor

    if sum_weights <= 0.0:
        fatal(
            f"Variable '{variable_name}', point index {point_index}, has non-positive "
            f"total weight."
        )

    #endif

    weighted_x = weighted_x_sum / sum_weights
    weighted_y = weighted_y_sum / sum_weights
    weighted_dy = 1.0 / math.sqrt(sum_weights)

    return weighted_x, weighted_y, weighted_dy


def print_progress_block(variable_name, point_index, input_paths, triples, combined_triple):
    """
    Print the three input period points and the resulting weighted average.
    """

    print("")
    print("------------------------------------------------------------")
    print(f"Variable: {variable_name}")
    print(f"Point index: {point_index}")
    print("")

    for input_index, triple in enumerate(triples):
        x_val, y_val, dy_val = triple

        print(
            f"  input {input_index + 1}: "
            f"{input_paths[input_index]}  "
            f"{{mean -t' = {x_val:.9f}, A = {y_val:.9f}, sigma = {dy_val:.9f}}}"
        )

    #endfor

    combined_x, combined_y, combined_dy = combined_triple

    print("")
    print(
        f"  combined:  "
        f"{{mean -t' = {combined_x:.9f}, A = {combined_y:.9f}, sigma = {combined_dy:.9f}}}"
    )
    print("------------------------------------------------------------")


def combine_files(file_data, input_paths):
    combined = OrderedDict()

    for variable_name in file_data[0].keys():
        n_points = len(file_data[0][variable_name])
        combined_triples = []

        for point_index in range(n_points):
            triples_to_combine = []

            for data in file_data:
                triples_to_combine.append(data[variable_name][point_index])

            #endfor

            combined_triple = weighted_average_triplet(
                triples_to_combine,
                variable_name,
                point_index
            )

            print_progress_block(
                variable_name,
                point_index,
                input_paths,
                triples_to_combine,
                combined_triple
            )

            combined_triples.append(combined_triple)

        #endfor

        combined[variable_name] = combined_triples

    #endfor

    return combined


def format_assignment(variable_name, triples):
    pieces = []

    for x_val, y_val, dy_val in triples:
        pieces.append(f"{{{x_val:.9f}, {y_val:.9f}, {dy_val:.9f}}}")

    #endfor

    return f"{variable_name} = {{{', '.join(pieces)}}};"


def write_output(path, combined):
    try:
        with open(path, "w") as output_file:
            for variable_name, triples in combined.items():
                output_file.write(format_assignment(variable_name, triples))
                output_file.write("\n")

            #endfor

    except OSError as exc:
        fatal(f"Could not write output file '{path}': {exc}")

    #endtry


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Combine three Mathematica-style asymmetry fit text files into one "
            "weighted-average output file."
        )
    )

    parser.add_argument("input_file_1", help="First input text file.")
    parser.add_argument("input_file_2", help="Second input text file.")
    parser.add_argument("input_file_3", help="Third input text file.")
    parser.add_argument("output_file", help="Output text file.")

    args = parser.parse_args()

    input_paths = [
        args.input_file_1,
        args.input_file_2,
        args.input_file_3
    ]

    file_data = []

    for path in input_paths:
        file_data.append(parse_file(path))

    #endfor

    validate_same_structure(file_data, input_paths)

    combined = combine_files(file_data, input_paths)

    write_output(args.output_file, combined)

    print("")
    print(f"Wrote weighted averages to: {args.output_file}")


if __name__ == "__main__":
    main()