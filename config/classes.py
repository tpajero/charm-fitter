class Parameter:
    """Class to define the objects related to a scan parameter in a single place."""

    def __init__(self, identifier, title, range1d, parametrisation, range2d=None, degrees=False):
        assert parametrisation in ["theo", "pheno", "both"]
        self.id = identifier
        self.title = title
        self.parametrisation = parametrisation
        self.range1d = range1d
        self.range2d = range2d if range2d else range1d
        self.degrees = degrees
        if self.degrees:
            self.title += r"\,[^\circ]$"


class Plot2d:
    """Class to define the properties of 2D plots."""

    def __init__(self, pars, logo="r", legpos="l"):
        self.pars = pars
        self.logo = logo
        self.legpos = legpos


class Subcombination:
    """Class to define the PDFs and title of a subcombination of measurements."""

    def __init__(self, title, pdfs):
        self.title = title
        self.pdfs = pdfs
