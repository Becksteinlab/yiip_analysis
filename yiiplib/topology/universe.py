import MDAnalysis as mda
import numpy as np

from MDAnalysis.core import topologyattrs

class SecondaryStructure(topologyattrs.ResidueAttr):
    attrname = "structures"
    singular = "structure"
    dtype = object

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.array(['' for _ in range(nr)], dtype=object)


class Helix(topologyattrs.ResidueAttr):
    attrname = "helices"
    singular = "helix"
    dtype = object

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.array(['' for _ in range(nr)], dtype=object)


class Domain(topologyattrs.ResidueAttr):
    attrname = "domains"
    singular = "domain"
    dtype = object

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.array(['' for _ in range(nr)], dtype=object)

class YiiPUniverse(mda.Universe):
    tmd_helices = [(9, 37),
                   (38, 66),
                   (79, 109),
                   (115, 142),
                   (149, 174),
                   (178, 208)]

    tmd_helix_names = ["M1",
                       "M2",
                       "M3",
                       "M4",
                       "M5",
                       "M6"]

    ctd_helices = [(213, 226),
                   (258, 277)]

    ctd_helix_names = ['M7',
                       'M8']

    ctd_sheets =  [(232, 241),
                   (244, 252),
                   (280, 287)]

    ctd_sheet_names = ['S1',
                       'S2',
                       'S3']

    @staticmethod
    def _selection_string_from_ranges(ranges):
        ids = [str(i) + "-" + str(j) for i, j in ranges]
        return "resid " + ' '.join(ids)

    @staticmethod
    def ctd_selection_string():
        helices = YiiPUniverse.ctd_helices + YiiPUniverse.ctd_sheets
        return YiiPUniverse._selection_string_from_ranges(helices)

    @staticmethod
    def tmd_selection_string():
        helices = YiiPUniverse.tmd_helices
        return YiiPUniverse._selection_string_from_ranges(helices)

    @staticmethod
    def _get_domain(resid):
        for lower, upper in YiiPUniverse.ctd_helices:
            if resid >= lower or resid <= upper:
                return "ctd"

        for lower, upper in YiiPUniverse.tmd_helices:
            if resid >= lower or resid <= upper:
                return "tmd"

        for lower, upper in YiiPUniverse.ctd_sheets:
            if resid >= lower or resid <= upper:
                return "tmd"

        return "neither"

    @staticmethod
    def _expand(i, j):
        return [j for j in range(i, j + 1)]

    def __init__(self, topology_file: str, **kwargs):
        super().__init__(topology_file)
        self.format = topology_file.split('.')[-1]
        self._bind_secondary_structure()
        self._bind_helix_names()
        self._bind_domain()

    def _bind_secondary_structure(self):
        """Assign either loop, sheet or helix to each residue.
        """
        values = ["" for _ in range(len(self.residues))]
        self.add_TopologyAttr(SecondaryStructure(values))

        prot = self.select_atoms("protein")
        prot.residues.structures = "loop"
        helices = self.ctd_helices + self.tmd_helices

        targets = np.concatenate([YiiPUniverse._expand(i, j) for i, j in helices])
        for t in targets:
            prot.residues[prot.residues.resids == t].structures = "helix"

        targets = np.concatenate([YiiPUniverse._expand(i, j) for i, j in self.ctd_sheets])
        for t in targets:
            prot.residues[prot.residues.resids == t].structures = "sheet"

    def _bind_helix_names(self):
        """Assign helix names to each residue.
        """
        values = ["" for _ in range(len(self.residues))]
        self.add_TopologyAttr(Helix(values))

        prot = self.select_atoms("protein")

        helices = [zip(self.ctd_helices, self.ctd_helix_names),
                   zip(self.ctd_helices, self.tmd_helix_names),
                   zip(self.ctd_sheets, self.ctd_sheet_names)]

        for hdefs in helices:
            for hrange, name in hdefs:
                targets = YiiPUniverse._expand(hrange[0], hrange[1])
                targets = np.array(targets, dtype=object)
                for t in targets:
                    prot.residues[prot.residues.resids == t].helices = name

    def _bind_domain(self):
        """Assign the domain to each residue.
        """
        values = ["" for _ in range(len(self.residues))]
        self.add_TopologyAttr(Domain(values))

        protomer = self.select_atoms("protein")

        for name in self.ctd_helix_names:
            protomer.residues[protomer.residues.helices == name].domains = "ctd"

        for name in self.tmd_helix_names:
            protomer.residues[protomer.residues.helices == name].domains = "tmd"

        for name in self.ctd_sheet_names:
            protomer.residues[protomer.residues.helices == name].domains = "ctd"

    def domain_select(self, protomer, domain):
        """Select an AtomGroup from a given domain in a given protomer.
        Parameters
        ----------
        protomer: str
            Protomer to select the atom group from. This can be either "A" or "B"
        domain: str
            Domain name. This can be either "tmd", "ctd", or "neither".
        Returns
        -------
        AtomGroup
            Atom group from selection
        """
        ag = self.select_atoms("protein and segid " + protomer)
        ag = ag.atoms[ag.atoms.domains == domain]
        return ag

    @property
    def ctd_A(self):
        return self.domain_select("A", "ctd")

    @property
    def ctd_B(self):
        return self.domain_select("B", "ctd")

    @property
    def tmd_A(self):
        return self.domain_select("A", "tmd")

    @property
    def tmd_B(self):
        return self.domain_select("B", "tmd")

    @property
    def neither_A(self):
        return self.domain_select("A", "neither")

    @property
    def neither_B(self):
        return self.domain_select("B", "neither")