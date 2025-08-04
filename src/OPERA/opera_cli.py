"""
Calling OPERA CLI from Python.
"""

import logging
import os
import glob
import libOPERA_Py as OPERA
	
class OPERACLI():

	def __init__(self):
		self.smiles_file_location = "src/Comparison/LocalIO/opera_input.smi" 
		self.predictions_full_path = "src/Comparison/LocalIO/Thermout.csv" 

		self.opera = OPERA.initialize()  # initializes OPERA instance

	def smiles_file(self, smiles_list, cid_list):
		"""
		Creates smiles file of molecules to be predicted.
		Inputs: smiles - list of SMILES.
		Outputs: smiles tempfile name.
		"""
		# Use LocalIO directory for consistency with system
		file = open(self.smiles_file_location, "w")
		string = ""
		for smiles, cid in zip(smiles_list, cid_list):
			string += smiles + "\t" + str(cid) + "\n"
		string = string[:-1]  # trims trailing newline

		file.write(string)
		file.close()
		return file

	def remove_file(self):
		"""
		Removes all temp files from LocalIO directory.
		"""
		for file in glob.glob(os.path.join(os.path.dirname(self.smiles_file_location), "opera_*")):
			try:
				os.remove(file)
			except Exception:
				pass  # Ignore errors if file doesn't exist or can't be removed

	def run_opera(self, smiles_list, cid_list, requested_properties=None):
		"""
		Runs OPERA CLI routine with dynamic property selection.
		
		Args:
			smiles_list: List of SMILES strings to process
			requested_properties: List of property keys to compute
		
		Returns:
			List of dictionaries containing only the requested property predictions
		"""
		try:
			self.remove_file()
			smiles_tempfile = self.smiles_file(smiles_list, cid_list)
			self.smiles_full_path = smiles_tempfile.name

			logging.warning("Initiating OPERA CLI execution.")
			args = ['-s', self.smiles_full_path, '-o', self.predictions_full_path, '-c'] + requested_properties
			self.opera.OPERA(*args)
			logging.warning("Finished executing OPERA CLI.")
			return True
			
		except Exception as e:
			logging.warning("Exception running OPERA: {}".format(e))
			self.remove_file()
			return False


