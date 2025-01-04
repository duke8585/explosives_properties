from oxygen_balance_emp import parse_sum_formula


# def detonation_products(molecule):
#     # Extract atom counts
#     C = molecule.get('C', 0)
#     H = molecule.get('H', 0)
#     N = molecule.get('N', 0)
#     O = molecule.get('O', 0)

#     # Initialize products
#     products = {'H2O': 0, 'CO2': 0, 'CO': 0, 'N2': 0, 'C(s)': 0}

#     # Track available oxygen
#     available_oxygen = O

#     # 1. Form H2O (2 H + 1 O)
#     h2o_count = min(H / 2, available_oxygen)
#     products['H2O'] = h2o_count
#     available_oxygen -= h2o_count
#     H -= h2o_count * 2  # Use up hydrogen

#     # 2. Form CO (1 C + 1 O)
#     co_count = min(C, available_oxygen)
#     products['CO'] = co_count
#     available_oxygen -= co_count
#     C -= co_count  # Use up carbon

#     # 3. Convert CO to CO2 if oxygen is still available (1 CO + 1 O → 1 CO2)
#     co2_count = min(products['CO'], available_oxygen)
#     products['CO2'] += co2_count  # Add to existing CO2 count
#     products['CO'] -= co2_count  # Remove converted CO
#     available_oxygen -= co2_count  # Use up oxygen

#     # 4. Form N2 (2 N atoms per molecule)
#     n2_count = N / 2  # Allow for fractional N2 molecules
#     products['N2'] = n2_count
#     N -= n2_count * 2  # Use up nitrogen

#     # 5. Remaining carbon becomes solid (C)
#     products['C(s)'] = C

#     return products


def detonation_products(atom_counts):
    """
    mechanism following following reaction hierarchy:
    1. Carbon atoms react with the available oxygen to form CO
    2. N atoms combine to form N2
    3. 1/3 of the CO produced is converted to C + CO2
    4. 1/6 of the original CO produced reacts with any available H in the compound to form C + H2O
    5. if there is still oxygen left, C and CO are converted to CO2
    see https://edu.rsc.org/uk-chemistry-olympiad/writing-explosive-equations-chemistry-olympiad-worked-answers/1062.article#:~:text=Thus%2C%20the%20oxygen%20balance%20of%20RDX%20is%20%E2%80%9321.6%25.
    """
    # Initialize counts for each element
    elements = {
        "C": atom_counts.get("C", 0),
        "O": atom_counts.get("O", 0),
        "N": atom_counts.get("N", 0),
        "H": atom_counts.get("H", 0),
    }

    print(f"Initial elements: {elements}")

    # Step 1: Carbon atoms react with the available oxygen to form CO
    CO = min(elements["C"], elements["O"])
    elements["C"] -= CO
    elements["O"] -= CO

    print(f"After forming CO: CO={CO:.3f}, elements={elements}")

    # Step 2: N atoms combine to form N2
    N2 = elements["N"] / 2
    elements["N"] -= N2 * 2

    print(f"After forming N2: N2={N2:.3f}, elements={elements}")

    # Step 3: 1/3 of the CO produced is converted to C + CO2
    CO2 = CO / 3
    C_from_CO2 = CO2
    CO -= CO2

    print(
        f"After converting CO to CO2: CO2={CO2:.3f}, CO={CO:.3f}, elements={elements}"
    )

    # Step 4: 1/6 of the original CO produced reacts with any available H in the compound to form C + H2O
    H2O = min(CO / 6, elements["H"] / 2)
    C_from_H2O = H2O
    CO -= H2O
    elements["H"] -= H2O * 2

    print(f"After forming H2O: H2O={H2O:.3f}, CO={CO:.3f}, elements={elements}")

    print(f"Remaining elements before STEP 5: {elements}")

    # Step 5: if there is still oxygen left, C and CO are converted to CO2
    if elements["O"] > 0:
        # Convert C_from_CO2 and C_from_H2O to CO2
        C_to_CO2 = min(C_from_CO2 + C_from_H2O, elements["O"] / 2)
        elements["O"] -= C_to_CO2 * 2
        CO2 += C_to_CO2

        print(
            f"After converting C_from_CO2 and C_from_H2O to CO2: CO2={CO2:.3f}, elements={elements}"
        )

        # Convert remaining C to CO2
        C_to_CO2 = min(elements["C"], elements["O"] / 2)
        elements["C"] -= C_to_CO2
        elements["O"] -= C_to_CO2 * 2
        CO2 += C_to_CO2

        print(
            f"After converting remaining C to CO2: CO2={CO2:.3f}, elements={elements}"
        )

        # Convert remaining CO to CO2
        CO_to_CO2 = min(CO, elements["O"])
        CO -= CO_to_CO2
        elements["O"] -= CO_to_CO2
        CO2 += CO_to_CO2

        print(
            f"After converting remaining CO to CO2: CO2={CO2:.3f}, CO={CO:.3f}, elements={elements}"
        )

    # Collect the detonation products
    products = {
        "CO": CO,
        "N2": N2,
        "CO2": CO2,
        "H2O": H2O,
        "C(s)": elements["C"],
        "O2": elements["O"] / 2,  # Remaining oxygen as O2
    }

    print(f"Final products: {products}")
    print("=" * 10, "\n" * 2)
    return products


#####
#####


def combust_products(molecule):
    # Extract atom counts
    C = molecule.get("C", 0)
    H = molecule.get("H", 0)
    N = molecule.get("N", 0)
    O = molecule.get("O", 0)

    # Initialize products
    products = {"H2O": 0, "CO2": 0, "N2": 0, "O2_excess": 0}

    # 1. Form H2O (2 H + 1 O)
    h2o_count = H / 2  # Fully combust all hydrogen
    products["H2O"] = h2o_count

    # 2. Form CO2 (1 C + 2 O)
    co2_count = C  # Fully combust all carbon
    products["CO2"] = co2_count

    # 3. Form N2 (2 N)
    n2_count = N / 2  # Pair up all nitrogen atoms
    products["N2"] = n2_count

    # 4. Calculate required oxygen (for H2O and CO2)
    required_oxygen = h2o_count + (co2_count * 2)

    # 5. Calculate excess oxygen (total oxygen in molecule minus used oxygen)
    available_oxygen = O
    excess_oxygen = (
        available_oxygen - required_oxygen if available_oxygen > required_oxygen else 0
    )

    # Include any excess oxygen or oxygen shortfall
    products["O2_excess"] = (
        required_oxygen - available_oxygen if available_oxygen < required_oxygen else 0
    )

    return products


if __name__ == "__main__":

    molecule = parse_sum_formula("C7H5N3O6")

    # Example molecule
    products = detonation_products(molecule)

    # Output the products
    print("Detonation products:")
    for product, count in products.items():
        if count > 0:
            print(f"{product}: {count:.2f}")

    # Example molecule
    products = combust_products(molecule)

    # Output the products
    print("Combustion products:")
    for product, count in products.items():
        if count > 0:
            print(f"{product}: {count:.2f}")
