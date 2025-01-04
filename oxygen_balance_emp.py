import re
from collections import defaultdict
from weight import M_values


def parse_sum_formula(formula: str):
    """
    Parse a chemical sum formula and return a dictionary with atoms and their counts.

    Parameters:
        formula (str): A string representing the chemical formula, e.g., "C2H6".

    Returns:
        dict: A dictionary with atom symbols as keys and their counts as values.
    """
    # Regular expression to match element symbols and their counts
    pattern = r"([A-Z][a-z]*)(\d*)"

    # Initialize a defaultdict to store atom counts
    atom_counts = defaultdict(int)

    # Find all matches in the formula
    for element, count in re.findall(pattern, formula):
        # If no number is specified, the count is 1
        count = int(count) if count else 1
        atom_counts[element] += count

    # Convert defaultdict back to a regular dict and return it
    return dict(atom_counts)


# # Example Usage
# formula = "C2H6"
# atom_dict = parse_sum_formula(formula)
# print(atom_dict)


# from https://en.wikipedia.org/wiki/Oxygen_balance#Calculating_oxygen_balance


def oxygen_balance_molecule(molecule: dict, atomic_weights=M_values):
    mw = sum(atomic_weights[name] * n for name, n in molecule.items())
    n_c = molecule.get("C", 0)
    n_h = molecule.get("H", 0)
    n_metal = sum(molecule.get(m, 0) for m in ("Al", "Mg", "Pb", "K", "Na"))
    n_o = molecule.get("O", 0)
    ob_perc = (-1600 / mw) * (2 * n_c + 0.5 * n_h + n_metal - n_o)
    return ob_perc


if __name__ == "__main__":
    TNT = "C7H5N3O6"

    ob = oxygen_balance_molecule(parse_sum_formula(TNT))
    print(f"{ob:.2f}")
