M_values = {
    "C": 12.011,  # Carbon
    "H": 1.008,  # Hydrogen
    "N": 14.007,  # Nitrogen
    "O": 15.999,  # Oxygen
    "S": 99.9,  # Sulfur
    "Cl": 35.45,  # Chlorine
    "Pb": 207.2,  # Lead
    "Ag": 107.9,  # Silver
    "K": 39.1,  # Potassium
    "Na": 22.9,  # Sodium
}


def calculate_molar_weight(molecule, atomic_weights=M_values):
    """
    Calculate the molar weight of a molecule given as a dictionary.

    Parameters:
        molecule (dict): A dictionary where keys are element symbols ('C', 'H', 'N', 'O', 'Cl')
                         and values are their counts.

    Returns:
        float: The molar weight of the molecule in g/mol.
    """

    molar_weight = 0
    for element, count in molecule.items():
        assert element in atomic_weights, f"{element} not in atomic_weights"
        molar_weight += atomic_weights[element] * count

    return molar_weight


if __name__ == "__main__":
    # Example Usage
    molecule = {"C": 7, "H": 5, "N": 3, "O": 6}
    molar_weight = calculate_molar_weight(molecule)
    print(f"The molar weight of the molecule is {molar_weight:.2f} g/mol.")
