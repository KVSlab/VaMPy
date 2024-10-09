#!/usr/bin/env python

# Program:   AneuTool
# Module:    NetworkBoundaryConditions.py
# Language:  Python
# Date:      $Date: 2016/17/04 00:00:00 $
# Version:   $Revision: 0.0.1 $
# Author:    Christophe Chnafa

#   Copyright (c) Christophe Chnafa. All rights reserved.


class FlowSplitting(object):
    """This class computes the flow splitting for a network according to
    [ADD REF AJRN] and [VALIDATION PAPER]"""

    def __init__(self):
        self.hasComputedAlphas = False
        self.hasComputedBetas = False
        self.hasComputedGammas = False

    def ComputeAlphas(self, network, verboseprint):
        """Compute the segments alpha coefficient.

        Compute for each blanked segment of a network its attributed alpha
        coefficient. For a network division with N daughter segments
        (N can be superior to two), the alpha coefficient of each segment i
        is computed as:
        alpha_i = S / sumSurfaces = S_i / sum_1^N(segment area_j)
        That means that alpha is the percentage of flow from the datum
        (the supplier) segment going though the i-th segment of the considered
        division. Segments that are not blanked are left with an alpha
        coefficient of 1.0 since there is no flow division.

        """
        if network.GetNumberOfOutlet() < 2:
            raise RuntimeError("The network has only one outlet.")
        # For each element of the network,
        for element in network.elements:
            if not (element.IsBlanked()):
                continue
            # find the adjacent blanked branches.
            adjacentBranchFound = False
            S = network.elements[element.GetFrontSegment()].GetMeanArea()
            sumSurfaces = S
            for otherElement in network.elements:
                if otherElement.GetId() == element.GetId():
                    continue
                if not (otherElement.IsBlanked()):
                    continue
                if otherElement.GetInPointsx0Id() == element.GetInPointsx0Id():
                    adjacentBranchFound = True
                    sumSurfaces += network.elements[
                        otherElement.GetFrontSegment()
                    ].GetMeanArea()

            # At least one adjacent blanked branch should be found.
            if not (adjacentBranchFound):
                raise RuntimeError(
                    "Unexpected error, "
                    "adjacent branch not found. Check the connectivity."
                )
            alpha = S / sumSurfaces
            element.SetAlpha(alpha)
        self.hasComputedAlphas = True

    def ComputeBetas(self, network, verboseprint):
        """Compute the outlets beta coefficient.

        Compute for each outlet segment of a network its attributed beta
        coefficient. This coefficient correspond to the percentage of
        flow going out of this segment from the inlet flow. It is the actual
        flow division at the outlets. Segments that are not outlets are left
        with a beta coefficient of 0.0. The algorithm compute the beta
        coefficients as the multiplications of the alpha coefficient from the
        considered outlet to the root of the network where the inlet is
        applied.

        """
        if not (self.hasComputedAlphas):
            raise RuntimeError("Alpha coefficients need to be computed first.")
        # For each element of the network, ...
        for element in network.elements:
            if not (element.IsAnOutlet()):
                continue
            foundNetworkRoot = False
            beta = 1.0
            currentElement = element
            # ... gather alpha coefficients at each division by scanning
            # the network in direction of the flow inlet.
            while not (foundNetworkRoot):
                if currentElement.GetBehindSegment() is None:
                    raise RuntimeError(
                        "The network is constitued of one segment "
                        "or the input centerlines have a hanging segment."
                    )
                behindElement = network.elements[currentElement.GetBehindSegment()]
                if behindElement.IsBlanked():
                    beta *= behindElement.GetAlpha()
                if behindElement.IsAnInlet():
                    foundNetworkRoot = True
                else:
                    currentElement = behindElement
            element.SetBeta(beta)
        self.hasComputedBetas = True

    def ComputeGammas(self, network, verboseprint):
        """Compute the outlets gamma coefficient.

        Compute for each outlet segment of a network its attributed beta
        coefficient. This coefficient correspond to the percentage of
        flow going out of this segment from the inlet flow. The Gamma
        coefficient of each outlet i is computed as:
        gamma_i = S_i / sumOutletAreas

        """
        sumAreas = 0.0
        for element in network.elements:
            if not (element.IsAnOutlet()):
                continue
            sumAreas += element.GetMeanArea()
        for element in network.elements:
            if not (element.IsAnOutlet()):
                continue
            gamma = element.GetMeanArea() / sumAreas
            element.SetGamma(gamma)
        self.hasComputedGammas = True

    def CheckTotalFlowRate(self, network, verboseprint):
        """Check if the sum of the outflows is 100%."""
        tol = 0.000001  # Flow balance error tolerance.
        if self.hasComputedBetas:
            sumBeta = 0.0
            for element in network.elements:
                if element.IsAnOutlet():
                    sumBeta += element.GetBeta()
            if abs(sumBeta - 1.0) > tol:
                raise RuntimeError("Unexpected error, sum(Beta) coefficients != 1.0")
        if self.hasComputedGammas:
            sumGamma = 0.0
            for element in network.elements:
                if element.IsAnOutlet():
                    sumGamma += element.GetGamma()
            if abs(sumGamma - 1.0) > tol:
                raise RuntimeError("Unexpected error, sum(Gamma) coefficients != 1.0")
