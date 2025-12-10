from tools.analyze import Analyze
from tools.calculation_functions import DifferenceFunction
import constants as c
import pandas as pd

pd.set_option('display.float_format', '{:.2e}'.format)

folder   = "/home/deck/Desktop/Bak_Baig/Program/Data"
solutes  = list(c.Filename_to_Solute.values())
solvents = list(c.Solvent_to_shorten.keys())
solutes_shorten  = list(c.Solute_to_shorten.values())
solvents_shorten = list(c.Solvent_to_shorten.values())

solvent_filename = "Differences_between_solvents.tex"
solute_filename  = "Differences_between_solutes.tex"

function = DifferenceFunction.symmetric_percentage_difference

analasys = Analyze(folder, function, solutes, solvents, solutes_shorten, solvents_shorten)
analasys.solute_differences(bdp_central="H", mapping=c.Atom_number_conversion)

# analasys.solvent_differences()
# # analasys.display_solvent_differences()
# # analasys.display_solvent_differences_by_solute()
# analasys.generate_latex_results_document(differences="Solvent",
#                                          use_solvent_by_solute=False)

# analasys.solute_differences(bdp_central="A")
# analasys.display_solute_differences()
# analasys.generate_latex_results_document(solute_filename, "Solute")

