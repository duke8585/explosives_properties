R = 8.31  # J/(mol*K)
p_std = 1.013e5  # Pa or N/m^2


def sum_gas_moles(molecules: dict):
    gases = ("CO2", "CO", "H2O", "N2", "O2", "H2")
    return sum(n for name, n in molecules.items() if name in gases)


def gas_volume(n: float = 1, T: float = 273):
    V = (n * R * T) / p_std  # V in m^3
    return V * 1000  # in l


from oxygen_balance_emp import parse_sum_formula
from combusion_products import detonation_products, combust_products
from weight import calculate_molar_weight
from oxygen_balance_emp import oxygen_balance_molecule


def detonation_volume_per_kg(molecule: dict, t: float):
    return gas_volume(sum_gas_moles(detonation_products(molecule)), t) / (
        calculate_molar_weight(molecule) * 1e-3
    )


def combustion_volume_per_kg(molecule: dict, t: float):
    return gas_volume(sum_gas_moles(combust_products(molecule)), t) / (
        calculate_molar_weight(molecule) * 1e-3
    )


explosives = {
    # nitroaromatics
    "TNT": "C7H5N3O6",
    "PA": "C6H3N3O7",
    "Tetryl": "C7H5N5O8",
    "HNS": "C14H6N6O12",
    # nitramines
    "RDX": "C3H6N6O6",
    "HMX": "C4H8N8O8",
    "CL20": "C6H6N12O12",
    # nitrate esters
    "EGDN": "C2H4N2O6",
    "Trinitroglycerin": "C3H5N3O9",
    "PETN": "C5H8N4O12",
    # nitroguanidines
    "NQ": "CH4N4O2",
    # ANFO
    "ANFO_3:1": "CN6H14O9",  # 10x for everything but C, because 95% AN, 5% F mixture, https://www.e-education.psu.edu/mng230/node/685 3 moles of nh4no3, 1 mole of fuel oil
    # others
    "Ammonium perchlorate": "NH4Cl1O4",
    "Ammonium chlorate": "NH4ClO3",
    "Tetrazene": "C2H8N10O",  # adjusted for *h2o
    "Pb(N3)2": "PbN6",
    "Ag3N": "Ag3N",
    "TATP": "C9H18O6",
    "NC_mono": "C6H9NO7",
    "NC_di": "C6H8N2O9",
    "NC_tri": "C6H7N3O11",
    "black powder": "K2N2O6SC3",  # 75% KNO3, 15% C, 10% S, https://www.quora.com/The-mass-ratio-of-black-powder-is-75-KNO3-15-C-10-S-If-KNO3-is-replaced-by-KClO3-what-is-the-ratio-of-black-powder
    "nitromethane": "CH3NO2",
    # reference
    "C(s)": "C",
    "C(s,poly)": "C100",
    "H2": "H2",
}

import csv

with open("explosives.csv", "w", newline="") as csvfile:
    fieldnames = [
        "name",
        "sum_formula",
        "M [g/mol]",
        "V_det_273K [l/kg]",
        "V_comb_273K [l/kg]",
        "OB [%]",
    ]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()

    for name, sum_formula in explosives.items():
        molecule = parse_sum_formula(sum_formula)
        print(f"processing {name}")
        mw = calculate_molar_weight(molecule)
        # print(
        #     f"detonation_volume_per_kg, 273K: {detonation_volume_per_kg(sum_formula, 273):.2f}"
        # )
        # print(
        #     f"combustion_volume_per_kg, 273K: {combustion_volume_per_kg(sum_formula, 273):.2f}"
        # )

        writer.writerow(
            {
                "name": name,
                "sum_formula": sum_formula,
                "M [g/mol]": mw,
                "V_det_273K [l/kg]": detonation_volume_per_kg(molecule, 273),
                "V_comb_273K [l/kg]": combustion_volume_per_kg(molecule, 273),
                "OB [%]": oxygen_balance_molecule(molecule),
            }
        )
