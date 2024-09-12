# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause


def normalizeValues(values):
    """
    Normalize a list of values to the range [0.45, 0.95], handling None values appropriately.

    Inputs:
    - values: A list of numerical values, which may include None values.

    Returns:
    - normalized_values: A list of normalized values, with None values preserved.
    """
    # Filter out None values from the input list
    non_none_values = [val for val in values if val is not None]

    # If there are no non-None values, return a list of None values of the same length
    if not non_none_values:
        return [None] * len(values)

    # Find the minimum and maximum values among the non-None values
    min_val = min(non_none_values)
    max_val = max(non_none_values)

    # If all non-None values are the same, return a list of 0.75 values
    if max_val == min_val:
        return [0.75] * len(values)

    # Normalize the values to the range [0.45, 0.95]
    normalized_values = [
        (
            None  # Keep None values as None
            if value is None
            else (0.5 + (value - min_val) / (max_val - min_val) * 0.5)
            - 0.05  # Normalize non-None values
        )
        for value in values
    ]

    return normalized_values


def normalizeDict(dictionary):
    """
    Normalize the values in a dictionary, where each key maps to a list of values.

    Inputs:
    - dictionary: A dictionary where keys are strings and values are lists of numerical values.

    Returns:
    - normalized_dict: A dictionary with the same structure as the input, but with normalized values.
    """
    # Determine the number of categories (length of the value lists)
    num_categories = len(dictionary[list(dictionary.keys())[0]])
    normalized_dict = {}

    # Iterate over each category index
    for category_index in range(num_categories):
        # Extract values for the current category from all dictionary entries
        category_values = [values[category_index] for values in dictionary.values()]

        # Normalize the extracted category values
        normalized_category = normalizeValues(category_values)

        # Assign the normalized values back to the corresponding keys in the dictionary
        for key, values in dictionary.items():
            if key not in normalized_dict:
                normalized_dict[key] = []
            normalized_dict[key].append(normalized_category.pop(0))

    return normalized_dict
