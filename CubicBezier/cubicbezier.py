import numpy as np
from numpy.polynomial import Polynomial as P
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from math import sqrt
from random import uniform as uni

class CubicBezier:
    """Class to represent a cubic bezier curve"""
    b0, b1, b2, b3 = P([1, -3, 3, -1]), P([0, 3, -6, 3]), P([0, 0, 3, -3]), P([0, 0, 0, 1])
    vb0, vb1, vb2, vb3 = b0.deriv(), b1.deriv(), b2.deriv(), b3.deriv()
    ab0, ab1, ab2, ab3 = vb0.deriv(), vb1.deriv(), vb2.deriv(), vb3.deriv()
    def __init__(self, x, y):
        """Initializes a cubic bezier curve

        Args:
            x ([float]): four y values
            y ([float]): four y values
        """
        self._nodes = np.array([x, y], np.float64).transpose()
    def velocity(self, t):
        """Velocity function of curve

        Args:
            t (float): 0 <= t <= 1

        Returns:
            numpy.ndarray: Normalized velocity vector
        """
        v = np.array([0, 0], np.float64)
        for i, p in zip(self._nodes, [CubicBezier.vb0, CubicBezier.vb1, CubicBezier.vb2, CubicBezier.vb3]):
            v += i * p(t)
        return v / np.linalg.norm(v)  # Normalize vector
    def position(self, t):
        """Position function of the curve

        Args:
            t (float): 0 <= t <= 1

        Returns:
            numpy.ndarray: Position vector
        """
        s = np.array([0, 0], np.float64)
        for i, p in zip(self._nodes, [CubicBezier.b0, CubicBezier.b1, CubicBezier.b2, CubicBezier.b3]):
            s += i * p(t)
        return s
    def bounding_box(self):
        """Calculates bounding box's corners of the curve

        Returns:
            [(float, float)]: List of four corners
                Bottom left, left top, right top, right bottom, ignored
        """
        a = -3*self._nodes[0] + 9*self._nodes[1] - 9*self._nodes[2] + 3*self._nodes[3]
        b =  6*self._nodes[0] -12*self._nodes[1] + 6*self._nodes[2]
        c = -3*self._nodes[0] + 3*self._nodes[1]

        dicsx = b[0] * b[0] - 4 * a[0] * c[0]
        dicsy = b[1] * b[1] - 4 * a[1] * c[1]
        x, y = [0, 1], [0, 1]
        if (dicsx >= 0):
            r1 = (-b[0] + sqrt(dicsx)) / (2*a[0])
            r2 = (-b[0] - sqrt(dicsx)) / (2*a[0])
            if r1 > 0 and r1 < 1:
                x.append(r1)
            if r2 > 0 and r2 < 1:
                x.append(r2)

        if (dicsy >= 0):
            r1 = (-b[1] + sqrt(dicsy)) / (2*a[1])
            r2 = (-b[1] - sqrt(dicsy)) / (2*a[1])
            if r1 > 0 and r1 < 1:
                y.append(r1)
            if r2 > 0 and r2 < 1:
                y.append(r2)

        maxx = max(self.position(i)[0] for i in x)
        minx = min(self.position(i)[0] for i in x)
        maxy = max(self.position(i)[1] for i in y)
        miny = min(self.position(i)[1] for i in y)

        return [(minx, miny), (minx, maxy), (maxx, maxy), (maxx, miny), (-1, -1)]

            
    def plot(self, xlim=(-0.1, 1.1), ylim=(-0.1, 1.1), nodes=False, tangents_and_normals=False, bounding_box=False):
        """Plots the curve with given optionals

        Args:
            xlim (tuple, optional): xlim of plot. Defaults to (-0.1, 1.1).
            ylim (tuple, optional): ylim of plot. Defaults to (-0.1, 1.1).
            nodes (bool, optional): Plot nodes that define curve. Defaults to False.
            tangents_and_normals (bool, optional): Plot tangent and normal vectors at equally spaced out t intervals. Defaults to False.
            bounding_box (bool, optional): Plot bounding box. Defaults to False.
        """
        
        # Set up plots
        Path = mpath.Path
        fig, ax = plt.subplots()
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], ylim[1])

        # Plot Curve
        pp1 = mpatches.PathPatch(
            Path(self._nodes,
                 [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
            fc="none", transform=ax.transData)
        ax.add_patch(pp1)

        if nodes:  # Plot nodes that define the curve
            for x, y in self._nodes:
                plt.plot(x, y, "go")

        if tangents_and_normals:  # Plot tangent and normal vectors
            for t in np.linspace(0.0, 1.0, num=10):
                v = self.velocity(t)
                origin = self.position(t)
                plt.quiver(*origin, v[0], v[1], color="r")
                plt.quiver(*origin, -v[1], v[0], color="b")

        if bounding_box:  # Plot bounding box
            pp2 = mpatches.PathPatch(
                Path(
                    self.bounding_box(),
                    [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]),
                    fc="none",
                    color="orange"
            )
            ax.add_patch(pp2)

        plt.show()

if __name__ == "__main__":
    x = CubicBezier([uni(0, 1), uni(0, 1), uni(0, 1), uni(0, 1)], 
                    [uni(0, 1), uni(0, 1), uni(0, 1), uni(0, 1)])
    x.plot(nodes=True, bounding_box=True)
    print(type(x.position(0.5)))