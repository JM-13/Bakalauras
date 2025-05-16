import re
import pandas as pd

#A massive speed up would be to save lines with program identifiers while reading the file for the first time and use that for each section. I am not writting that
#A small speed up would be to create lists of all the identifiers on class creation. I am not writting that
#A more cleaver way to return RAD lists without line numbers would slightly speed it up. I am not writting that
#The program will break if RAD is retrieved without line numbers and then Scan retrieval is run. I am not fixing that
class Extract:
    def __init__(self, filename=""):
        self.file_programs = {}

        self.filename = filename
        self.File_Lines = [(0, ""),]

        if self.filename:
            with open(self.filename, 'r') as infile:
                for line_number, line in enumerate(infile, start=1):
                    self.File_Lines.append((line_number, line))
                    if " (Enter " in line:
                        program = line.rstrip(')').split('/')[-1].split('.')[0]
                        if program not in list(self.file_programs.keys()):
                            self.file_programs[program] = []
                        self.file_programs[program].append(line_number)

        self.File_len = len(self.File_Lines)

        self.Energies = {}
        self.Energies_S0 = {}
        self.Energies_S1 = {}

        self.Coordinates = []
        self.Coordinates_XYZ = []
        self.Coordinates_RAD = []
        self.Coordinates_Optimized_RAD = []

        self.Scan_Fixed_Value = []
        self.Scan_Energies_S0 = []
        self.Scan_Energies_S1 = []
        self.Scan_Fixed_Coordinate_values = []
        self.Scan_Fixed_Coordinate_averages = []
        # self.Scan_Fixed_Coordinate_other = []

    def Energy_values(self, Program_identifier, Line_identifier): #Base function for energy value retrieval (not intended to be run on its own)
        #Program_identifier - what the program in the file is called
        #Line_identifier - text that indentifies the correct line
        self.Energies = {}

        Match_E_value = re.compile(r"[-]?\d+\.\d+") #Match a float in a string: [-]? optional minus, \d+ one or more didgets, \. period, \d+ one or more didgets after period

        for program_start in self.file_programs[Program_identifier]:
            for i in range(program_start, self.File_len):
                line_number, line = self.File_Lines[i]

                if Line_identifier in line: #If identifeing text in line
                        E_string_match = re.search(Match_E_value, line)
                        E_val = float(E_string_match.group()) #Extract energy value
                        self.Energies[line_number] = E_val
                        break

                if 'Leave Link' in line: #If program (in file text) ended. Fallback if file structure is not as anticipated
                        break

    def Coordinate_values(self, Program_identifier, Line_start_identifier, Line_end_identifier, Column_names, Row_value_types, Skip_lines = 0, Needed_row_elements = [], with_line_numbers=False): #Base function for coordinate retrieval (not intended to be run on its own)
        #Program_identifier       - what the program in the file is called
        #Line_start_identifier    - text that indentifies where the table starts
        #Line_end_identifier      - text that indentifies where the table ends
        #Column_names             - names for the file table columns (can be any names)
        #Row_value_types          - types of values that are in the row of the file table
        #Skip_lines = 0           - how many lines to pass after "Line_start_identifier" before starting to save them
        #Needed_row_elements = [] - which values in the row of the file table are needed; if it is set then "Column_names" and "Row_value_types" need to be set accordingly
        #with_line_numbers=False  - if the line numbers from which the rows are saved should also be saved

        self.Coordinates = []
        Table = [] #Stores rows before they are turned into a dataframe

        Row_element_number = len(Column_names) #How many values to expect
        if with_line_numbers:
            Column_names = ['Line Number'] + Column_names
        Skip = Skip_lines


        Save_as_row = False
        for program_start in self.file_programs[Program_identifier]:
            for i in range(program_start, self.File_len):
                line_number, line = self.File_Lines[i]

                if Save_as_row:
                    if Skip > 0:
                        Skip -= 1
                        continue
                    if Line_end_identifier in line:
                        df = pd.DataFrame(Table, columns=Column_names)
                        self.Coordinates.append(df)

                        Save_as_row = False
                        Skip = Skip_lines
                        Table = []
                        break

                    raw_row = line.split() #Convert string into a list; spaces in string differentiate between values
                    if Needed_row_elements:
                        row = [raw_row[i] for i in Needed_row_elements]
                    else:
                        row = raw_row

                    for elem in range(Row_element_number):
                        row[elem] = Row_value_types[elem](row[elem]) #Convert elements in list to the given types

                    if with_line_numbers:
                        row = [line_number] + row

                    Table.append(row)
                    continue

                if Line_start_identifier in line:
                    Save_as_row = True
                    continue

                if 'Leave Link' in line:
                    break

    def S0(self, with_line_numbers=False): #Get base energy S0
        Program_identifier = 'l502'
        Line_identifier = 'SCF Done:'

        self.Energy_values(Program_identifier, Line_identifier)
        self.Energies_S0 = self.Energies

        if with_line_numbers:
            return self.Energies_S0
        else:
            return list(self.Energies_S0.values())


    def S1(self, with_line_numbers=False): #Get excited energy S1
        Program_identifier = 'l914'
        Line_identifier = 'Total Energy'

        self.Energy_values(Program_identifier, Line_identifier)
        self.Energies_S1 = self.Energies

        if with_line_numbers:
            return self.Energies_S1
        else:
            return list(self.Energies_S1.values())

    def XYZ(self, with_line_numbers=False): #Get XYZ coordinate values of atoms
        Program_identifier = 'l202'
        Line_start_identifier = ' Number     Number       Type             X           Y           Z'
        Line_end_identifier = ' ---------------------------------------------------------------------'
        Column_names = ['Center Number', 'Atomic Number', 'Atomic Type', 'X', 'Y', 'Z']
        Row_value_types = [int, int, int, float, float, float]
        Skip_lines = 1

        self.Coordinate_values(Program_identifier, Line_start_identifier, Line_end_identifier, Column_names, Row_value_types, Skip_lines, with_line_numbers = with_line_numbers)
        self.Coordinates_XYZ = self.Coordinates

        return self.Coordinates_XYZ

    def RAD(self, with_line_numbers=False): #Get RAD (inner) coordinate values of atoms
        Program_identifier = 'l103'
        Line_start_identifier = '                                 (Linear)    (Quad)   (Total)'
        Line_end_identifier = 'Converged?'
        Column_names = ['Variable', 'Old X', '-DE/DX', 'Delta X (Linear)', 'Delta X (Quad)', 'Delta X (Total)', 'New X']
        Row_value_types = [str, float, float, float, float, float, float]

        self.Coordinate_values(Program_identifier, Line_start_identifier, Line_end_identifier, Column_names, Row_value_types, with_line_numbers = with_line_numbers)
        self.Coordinates_RAD = self.Coordinates

        return self.Coordinates_RAD

    def Optimized_RAD(self, with_line_numbers=False): #Get RAD (inner) coordinate values of atoms but from different table (it round one more, but has atom numbering)
        Program_identifier = 'l103'
        Line_start_identifier = '                           !   Optimized Parameters   !'
        Line_end_identifier = ' --------------------------------------------------------------------------------'
        Column_names = ['Name', 'Definition', 'Value', 'Derivative Info(-DE/DX)']
        Row_value_types = [str, str, float, float]
        Skip_lines = 4
        Needed_row_elements = [1, 2, 3, 6]

        self.Coordinate_values(Program_identifier, Line_start_identifier, Line_end_identifier, Column_names, Row_value_types, Skip_lines, Needed_row_elements, with_line_numbers)
        self.Coordinates_Optimized_RAD = self.Coordinates

        return self.Coordinates_Optimized_RAD

    def Scan_fixed(self, filename="", with_line_numbers=False): #Get the fixed coordinate of the Scan file
        Line_start_identifier = ' The following ModRedundant input section has been read:'

        Fixed_value = None
        Save = False

        if filename:
            with open(filename, 'r') as infile:
                for line_number, line in enumerate(infile, start=1):

                    if Line_start_identifier in line:
                        Save = True
                        continue

                    if Save:
                        Fixed_value = line.split()
                        Save = False
                        break
        else:
            for line_number, line in self.File_Lines:

                if Line_start_identifier in line:
                    Save = True
                    continue

                if Save:
                    Fixed_value = line.split()
                    Save = False
                    break

        Fixed_value_name = Fixed_value[0] + '(' + ",".join(Fixed_value[1:5]) + ')'

        self.Scan_Fixed_Value = [line_number, [Fixed_value_name] + Fixed_value[5:8]]

        if with_line_numbers:
            return self.Scan_Fixed_Value
        else:
            return self.Scan_Fixed_Value[1]

    def Scan_Optimized_Energy(self): #Get the S0 and S1 energys of every fixed coordinate incrament

        if not self.Energies_S0:
            S0 = self.S0(with_line_numbers = True)
        else:
            S0 = self.Energies_S0
        if not self.Energies_S1:
            S1 = self.S1(with_line_numbers = True)
        else:
            S1 = self.Energies_S1
        if not self.Coordinates_Optimized_RAD:
            RAD = self.Optimized_RAD(with_line_numbers = True)
        else:
            RAD = self.Coordinates_Optimized_RAD


        RAD_start_line_numbers = []
        for df in RAD:
            RAD_start_line_numbers.append(df.iloc[0, 0])

        S0_line_numbers = list(S0.keys())
        S1_line_numbers = list(S1.keys())
        for rad_l_n in RAD_start_line_numbers:

            S0_left_of_at = 0
            for s0_l_n in S0_line_numbers[S0_left_of_at:]:
                if rad_l_n>s0_l_n:
                    Optimized_S0 = s0_l_n
                    S0_left_of_at+=1
                else:
                    break
            self.Scan_Energies_S0.append(S0[Optimized_S0])

            S1_left_of_at = 0
            for s1_l_n in S1_line_numbers[S1_left_of_at:]:
                if rad_l_n>s1_l_n:
                    Optimized_S1 = s1_l_n
                    S1_left_of_at+=1
                else:
                    break
            self.Scan_Energies_S1.append(S1[Optimized_S1])

        return [self.Scan_Energies_S0, self.Scan_Energies_S1]

    def Scan_RAD_fixed_values(self): #Get the value of every fixed coordinate incrament

        if not self.Coordinates_Optimized_RAD:
            RAD = self.Optimized_RAD(with_line_numbers = True)
        else:
            RAD = self.Coordinates_Optimized_RAD
        if not self.Scan_Fixed_Value:
            Fixed_value = self.Scan_fixed()
        else:
            Fixed_value = self.Scan_Fixed_Value[1]


        Row_with_fixed_value = RAD[0].loc[RAD[0]['Definition'] == Fixed_value[0]].index[0] #Get the index of the table row that has the fixed coordinate
        for df in RAD:
            self.Scan_Fixed_Coordinate_values.append(df.iloc[Row_with_fixed_value]['Value'])

        return self.Scan_Fixed_Coordinate_values

    def Scan_RAD_fixed_average(self, return_fixed_opposit=False):

        if not self.Coordinates_Optimized_RAD: #Get the value of the average of every fixed coordinate incrament and its opposit
            RAD = self.Optimized_RAD(with_line_numbers = True)
        else:
            RAD = self.Coordinates_Optimized_RAD
        if not self.Scan_Fixed_Value:
            Fixed_value = self.Scan_fixed()
        else:
            Fixed_value = self.Scan_Fixed_Value[1]
        if not self.Scan_Fixed_Coordinate_values:
            RAD_fixed_values = self.Scan_RAD_fixed_values()
        else:
            RAD_fixed_values = self.Scan_Fixed_Coordinate_values


        Atom_listing_pattern = re.compile(r'D\((\d+),(\d+),(\d+),(\d+)\)')  #Match atoms in pattern 'D(a1,a2,a3,a4)', where a* are atom numbers
        a1, a2, a3, a4 = re.match(Atom_listing_pattern, Fixed_value[0]).groups()
        Atom_patterns = [re.compile(fr"D\(\d+,{a2},{a3},\d+\)"),
                         re.compile(fr"D\({a1},{a2},{a3},\d+\)"),
                         re.compile(fr"D\(\d+,{a2},{a3},{a4}\)")]

        RAD_fixed_values_opposit = []
        Fixed_value_opposit = []
        for df in RAD:
            df_filltered = df[df['Definition'].str.contains(Atom_patterns[0], regex=True)] #Which rows match this atom pattern; there will be 4
            df_filltered = df_filltered[~df_filltered['Definition'].str.contains(Atom_patterns[1], regex=True)] #Which rows do not match this atom pattern; will leave 2
            df_filltered = df_filltered[~df_filltered['Definition'].str.contains(Atom_patterns[2], regex=True)] #Which rows do not match this atom pattern; will leave the wanted row

            RAD_fixed_values_opposit.append(df_filltered.iloc[0]['Value']) #Extract opposit to fixed rows coordinate value

            if not Fixed_value_opposit:
                Fixed_value_opposit = df_filltered.iloc[0].tolist()[2] #Save the opposit to fixed pattern

        for fixd, opps in zip(RAD_fixed_values, RAD_fixed_values_opposit):
            self.Scan_Fixed_Coordinate_averages.append(round((fixd+opps)/2, 5)) #Save the average of the fixed and opposit to it value correctly? rounded

        if return_fixed_opposit:
            return [Fixed_value_opposit, self.Scan_Fixed_Coordinate_averages]
        else:
            return self.Scan_Fixed_Coordinate_averages

    # def Scan_RAD_other(self):
    #
    #     Other_coordinates = []
    #     for df in RAD:
    #         Other_coordinates.append(df.drop([Row_with_fixed_value, Row_with_fixed_value_opposit]))




















