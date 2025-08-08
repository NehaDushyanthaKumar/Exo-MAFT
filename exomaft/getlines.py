import os

class Linelist:
    """
    A class to process and filter atomic/molecular spectral line lists for specific elements or molecules,
    typically used in exoplanet transmission spectroscopy analysis.

    Attributes
    ----------
    elem_dict : dict
        Dictionary mapping element/molecule names to their atomic/molecular numbers for identification.
    elem : str
        The element/molecule to extract from the linelist.
    lines_dir : str
        Path to the base directory containing linelist files.
    outlines_dir : str
        Path to the directory containing raw line list text files.
    outlines_sorted_dir : str
        Path to the directory where cleaned and sorted files will be saved.
    outlines_master : str
        Path to the final merged line list file.
    """

    elem_dict = {
        "H20": 106.0,
        # You can extend this with real atomic/molecular numbers:
        # "CH4": ...,
        # "CO": ...,
        # "NH3": ...,
    }

    def __init__(self, elem):
        """
        Initialize the Linelist processor for a given element/molecule.

        Parameters
        ----------
        elem : str
            The key from `elem_dict` representing the target species (e.g., "H20").
        """
        self.elem = elem
        self.lines_dir = "exomaft/data"
        self.outlines_dir = os.path.join(self.lines_dir, "outlines/")
        self.outlines_sorted_dir = os.path.join(self.lines_dir, "outlines_sorted/")
        self.outlines_master = os.path.join(self.lines_dir, "outlines_master.txt")

        print(self.outlines_dir)

    def sortlines(self):
        """
        Filter, clean, and sort individual outline linelist files for the specified element.

        This method:
        - Iterates through all `.txt` files in `outlines_dir`.
        - Parses each line to extract relevant columns.
        - Filters lines based on atomic/molecular number matching the selected element.
        - Sorts the cleaned data by wavelength.
        - Writes cleaned lines to new `_cleaned.txt` files in `outlines_sorted_dir`.
        """
        for filename in os.listdir(self.outlines_dir):

            if not filename.endswith(".txt"):
                continue

            file_path = os.path.join(self.outlines_dir, filename)
            print(file_path)

            clean_rows = []

            with open(file_path, "r", encoding="ascii") as f:
                lines = f.readlines()

                for line in lines[1:]:  # skip header
                    columns = line.strip().split()

                    try:
                        col1 = abs(float(columns[0]))   # Wavelength
                        col2 = float(columns[1])        # Atomic/molecular number
                        col3 = columns[2]               # Energy (eV)
                        col4 = columns[3]               # log(gf)
                    except Exception as e:
                        print(f"ERROR parsing line: {repr(line)}")
                        print(f"Exception: {e}")
                        raise

                    try:
                        if str(col2).startswith(str(self.elem_dict[self.elem])):
                            print(f"{self.elem} line at {col1}")
                            clean_rows.append((col1, col2, col3, col4))
                    except:
                        print(f"No {self.elem} in wavelength range.")
                        return

                clean_rows.sort()

            cleaned_filename = filename.replace(".txt", "") + "_cleaned.txt"

            with open(os.path.join(self.outlines_sorted_dir, cleaned_filename), "w") as f:
                for col1, col2, col3, col4 in clean_rows:
                    f.write(f"{col1} {col2} {col3} {col4}\n")

            print(f"Wrote cleaned and sorted line list to {cleaned_filename}")

    def combine(self):
        """
        Combine all cleaned linelist files into a single master linelist file.

        This method:
        - Iterates through all cleaned `.txt` files in `outlines_sorted_dir`.
        - Merges them into one master file located at `outlines_master`.
        """
        with open(self.outlines_master, 'w') as outfile:
            for filename in sorted(os.listdir(self.outlines_sorted_dir)):
                if filename.endswith(".txt"):
                    filepath = os.path.join(self.outlines_sorted_dir, filename)
                    with open(filepath, 'r') as infile:
                        for line in infile:
                            outfile.write(line)

        print(f"Master outfile written to {self.outlines_master}")
