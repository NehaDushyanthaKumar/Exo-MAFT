import os 

class Linelist:

    elem_dict= {
        "H20": 106.0,
        #"CH4": ,
        #"CO": ,
        #"NH3": , 
    }

    def __init__(self,elem):
        self.elem = elem
        self.lines_dir = "/Users/evan/codeastro/group_project/data/"
        self.outlines_dir = os.path.join(self.lines_dir, "outlines/")
        self.outlines_sorted_dir = os.path.join(self.lines_dir, "outlines_sorted/")
        self.outlines_master = os.path.join(self.lines_dir, "outlines_master.txt")

        print(self.outlines_dir)

    def sortlines(self):

        #open outlines file and loop through every file inside
        for filename in os.listdir(self.outlines_dir):

            if not filename.endswith(".txt"):  # or your expected extension
                continue

            file_path = os.path.join(self.outlines_dir,filename)

            print(file_path)

            clean_rows = []

            #open line list file in outlines
            with open(file_path,"r",encoding="ascii") as f:

                lines = f.readlines()

                for line in lines[1:]:
                    #split columns [wavelength, atomic num, eV, loggf, damp const, flag] -- we only want the first 4
                    columns = line.strip().split()

                    try:
                        col1 = abs(float(columns[0]))
                        col2 = float(columns[1])
                        col3 = columns[2]
                        col4 = columns[3]                    
                    except Exception as e:
                        print(f"ERROR parsing line: {repr(line)}")
                        print(f"Exception: {e}")
                        raise

                
                    #only keep rows with the molecule or element we want 
                    try: 
                        if str(col2).startswith(str(self.elem_dict[self.elem])):
                            print(f"{self.elem} line at {col1}")
                            clean_rows.append((col1,col2,col3,col4))

                    except:
                        print(f"No {self.elem} in wavelength range.")
                        return

                clean_rows.sort()

            cleaned_filename = filename.replace(".txt", "") + "_cleaned.txt"

            with open(self.outlines_sorted_dir + cleaned_filename,"w") as f:
                for col1, col2, col3, col4 in clean_rows:
                    f.write(f"{col1} {col2} {col3} {col4}\n")
                print(f'Wrote cleaned and sorted line list to {filename.replace(".txt", "") + "_cleaned.txt"}')


    def combine(self):
        with open(self.outlines_master, 'w') as outfile:
            for filename in sorted(os.listdir(self.outlines_sorted_dir)):
                if filename.endswith(".txt"):
                    filepath = os.path.join(self.outlines_sorted_dir,filename)
                    with open(filepath,'r') as infile:
                        for line in infile:
                            outfile.write(line) 
        print(f"Master outfile written to {self.outlines_master}")
