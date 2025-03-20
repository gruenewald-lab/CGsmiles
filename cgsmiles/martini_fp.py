from .martini_fp_base import MartiniFingerPrinterBase

class MatrixFingerPrint(MartiniFingerPrinterBase):
    """
    Same as in CGsmiles paper.
    """

    def __init__(self):
        super().__init__()

    def aggregate_features(self, cgmol):
        fp_matrix = np.zeros((20, 6))
        for node in cgmol.nodes:
            size_arr = cgmol.nodes[node]['size']
            level_arr = cgmol.nodes[node]['inter_array']
            fp_matrix[level_arr, size_arr] += 1
        return fp_matrix

    def generate_fingerprint(self, cgmol, solvmol):
        targets = ["W"]
        targets += list(nx.get_node_attributes(solvmol, "fragname").values())
        self.annotate_features(self, cg_mol, targets, label_name='fragname')
        fp_matrix = self.aggregate_features(self, cgmol)
        return fp_matrix
