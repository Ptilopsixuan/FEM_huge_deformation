import numpy as np
import math

class Point: # 节点
    id = 0
    # coordinate
    x = 0
    y = 0
    # load
    px = 0
    py = 0
    pm = 0
    # displacement
    ax = 0
    ay = 0
    # axial displacement
    u = 0
    # deflection
    w = 0
    # rotation
    theta = 0

    def __init__(self, data:list[float]):
        (self.id, self.x, self.y) = data

    def __str__(self):
        return f"Node {self.id}: ({self.x}, {self.y})"
    
    def __repr__(self):
        return self.__str__()

class Material: # 材料
    id = 0
    E = 0
    v = 0
    G = 0

    def __init__(self, data: list[float]) -> None:
        (self.id, self.E, self.v) = data
        # calculate shear modulus
        self.G = self.E / (2 * (1 + self.v))

    def __str__(self) -> str:
        return f"<id:{self.id}, E:{self.E}>, G:{self.G}, v:{self.v}>"

    def __repr__(self) -> str:
        return self.__str__()

class Plane: # 截面
    id = 0
    b = 0
    h = 0
    # area
    A = 0
    # moment of inertia
    I = 0
    # shape of section
    '''circle = 0, rectangle = 1'''
    shape = 0

    def __init__(self,  data: list[float]) -> None:
        (self.id, self.shape, self.b, self.h) = data
        # calculate area
        self.A = self.b * self.h
        # calculate moment of inertia
        if self.shape == 0: # circle
            self.I = math.pi * self.b**4 / 64
        elif self.shape == 1: # rectangle
            self.I = self.b * self.h**3 / 12

    def __str__(self) -> str:
        return f"<id:{self.id}, A:{self.A}, I:{self.I}, shape:{self.shape}>"

    def __repr__(self) -> str:
        return self.__str__()

class Unit: # 单元
    #region parameters
    id = 0
    # node id
    i = 0
    j = 0
    # material
    material: Material = None
    # plane
    plane: Plane = None
    # length
    L = 0
    # angle
    theta = 0
    # trsnformation matrix
    c = 0
    s = 0
    T = np.zeros((6, 6))
    # axial stiffness
    EA = 0
    # flexural stiffness
    EI = 0
    # shear stiffness
    GA = 0
    # Correction factor for shear distribution
    k:float = 0
    b:float = 0
    # internal force
    i_N = 0
    i_V = 0
    i_M = 0
    j_N = 0
    j_V = 0
    j_M = 0
    # endregion

    def __init__(self, data:list[float], shear: bool = False) -> None:
        (self.id, self.i, self.j, self.material, self.plane) = data
        # calculate length
        dx = self.j.x - self.i.x
        dy = self.j.y - self.i.y
        self.L = np.hypot(dx, dy)
        # calculate angle
        dx = self.j.x - self.i.x
        dy = self.j.y - self.i.y
        self.theta = np.arctan2(dy, dx)
        # calculate axial stiffness
        self.EA = self.material.E * self.plane.A
        # calculate flexural stiffness
        self.EI = self.material.E * self.plane.I
        # calculate shear stiffness
        self.GA = self.material.G * self.plane.A
        # calculate stiffness matrix
        if shear:
            self.Ke = self.calculateKe_with_shear()
        else:
            self.Ke = self.calculateKe_without_shear()

    def __str__(self):
        return f"Element {self.id}: ({self.i}, {self.j})"

    def __repr__(self):
        return self.__str__()
    
    def calculateT(self) -> np.ndarray:
        # calculate transformation matrix
        c = np.cos(self.theta)
        s = np.sin(self.theta)
        self.T = np.array([[c, s, 0, 0, 0, 0],
                           [-s, c, 0, 0, 0, 0],
                           [0, 0, 1, 0, 0, 0],
                           [0, 0, 0, c, s, 0],
                           [0, 0, 0, -s, c, 0],
                           [0, 0, 0, 0, 0, 1]])
        return self.T
    
    def calculateKe_without_shear(self) -> np.ndarray:
        # without shear deformation
        k1 = self.EA / self.L
        k2 = 12 * self.EI / self.L**3
        k3 = 6 * self.EI / self.L**2
        k4 = 4 * self.EI / self.L
        k5 = 2 * self.EI / self.L
        Ke = np.array([[k1, 0, 0, -k1, 0, 0],
                       [0, k2, k3, 0, -k2, k3],
                       [0, k3, k4, 0, -k3, k5],
                       [-k1, 0, 0, k1, 0, 0],
                       [0, -k2, -k3, 0, k2, -k3],
                       [0, k3, k5, 0, -k3, k4]])
        return Ke
    
    def calculateKe_with_shear(self) -> np.ndarray:
        # with shear deformation
        self.k = [10/9 if self.plane.shape == 0 else 6/5] # P312
        self.b = 12 * self.EI * self.k / self.GA / self.L**2
        k0 = self.EI / (1 + self.b) / self.L**3
        k1 = self.EA / self.L / k0
        Ke = np.array([ [k1, 0, 0, -k1, 0, 0],
                        [0, 12, 6 * self.L, 0, -12, 6 * self.L],
                        [0, 6 * self.L, (4 + self.b) * self.L**2, 0, -6 * self.L, (2 - self.b) * self.L**2],
                        [-k1, 0, 0, k1, 0, 0],
                        [0, -12, -6 * self.L, 0, 12, -6 * self.L],
                        [0, 6 * self.L, (2 - self.b) * self.L**2, 0, -6 * self.L, (4 + self.b) * self.L**2]])
        Ke = k0 * Ke
        return Ke

    def transform(self) -> np.ndarray:
        # transform element stiffness matrix
        self.calculateT()
        self.Ke_ba = np.dot(np.dot(self.T.T, self.Ke), self.T)
        return self.Ke_ba

class Payload: # 荷载
    point: Point = None
    px = 0
    py = 0
    pm = 0

    def __init__(self, data: list) -> None:
        (self.point, self.px, self.py, self.pm) = data

    def __str__(self) -> str:
        return f"<point:{self.point.id}, px:{self.px}, py:{self.py}, pm:{self.pm}>"

    def __repr__(self) -> str:
        return self.__str__()

class Constraint: # 约束
    point: Point = None
    axis = 1 # 1: x, 2: y, 3: theta
    value = 0

    def __init__(self, data: list) -> None:
        (self.point, self.axis, self.value) = data

    def __str__(self) -> str:
        return f"<point:{self.point.id}, axis:{self.axis}, value:{self.value}>"

    def __repr__(self) -> str:
        return self.__str__()
    
class Hitch: # 铰接
    unit: Unit = None
    points: list[Point] = []

    def __init__(self, unit: Unit) -> None:
        self.unit = unit
    def __str__(self) -> str:
        return f"<unit:{self.unit.id}, point:{[p.id for p in self.points]}>"
    def __repr__(self) -> str:
        return self.__str__()
    def addPoint(self, point: Point) -> None:
        self.points.append(point)

class InputData:
    points: list[Point] = []
    planes: list[Plane] = []
    units: list[Unit] = []
    materials: list[Material] = []
    payloads: list[Payload] = []
    constraints: list[Constraint] = []
    hitches: list[Hitch] = []

    def __init__(self) -> None:
        self.points, self.planes, self.units, self.materials, self.payloads, self.constraints = [], [], [], [], [], []

    def __str__(self) -> str:
        return f"<\npoints:{self.points},\nplanes:{self.planes},\nunits:{self.units},\nmaterials:{self.materials},\npayloads:{self.payloads},\nconstraints:{self.constraints}\n>"

    def __repr__(self) -> str:
        return self.__str__()
    
    def addConstraint(self) -> list[int]:
        ConstraintLine = []
        for constraint in self.constraints:
            ConstraintLine.append(int((constraint.point.id-1)*3+constraint.axis-1))

        return ConstraintLine

    def calculateHitch(self, unit: Unit) -> None:
        for hitch in self.hitches:
            if hitch.unit.id == unit.id:
                index = self.units.index(unit)
                if len(hitch.points) == 1:
                    # replace_k = 3 * unit.EI / unit.L**3
                    if hitch.points[0].id == hitch.unit.i.id:
                        self.units[index].Ke[1][1] = self.units[index].Ke[4][4] = self.units[index].Ke[1][1] / 4
                        self.units[index].Ke[4][1] = self.units[index].Ke[1][4] = self.units[index].Ke[4][1] / 4
                        self.units[index].Ke[1][5] = self.units[index].Ke[5][1] = self.units[index].Ke[1][5] / 2
                        self.units[index].Ke[4][5] = self.units[index].Ke[5][4] = self.units[index].Ke[4][5] / 2
                        self.units[index].Ke[5][5] = (self.units[index].Ke[5][5] - self.units[index].b) * 3 / 4 + self.units[index].b
                        self.units[index].Ke[1][2] = self.units[index].Ke[2][4] = self.units[index].Ke[2][5] = self.units[index].Ke[2][2] = \
                        self.units[index].Ke[5][2] = self.units[index].Ke[4][2] = self.units[index].Ke[2][1] = 0
                    elif hitch.points[0].id == hitch.unit.j.id:
                        self.units[index].Ke[1][1] = self.units[index].Ke[4][4] = self.units[index].Ke[1][1] / 4
                        self.units[index].Ke[4][1] = self.units[index].Ke[1][4] = self.units[index].Ke[4][1] / 4
                        self.units[index].Ke[1][2] = self.units[index].Ke[2][1] = self.units[index].Ke[1][2] / 2
                        self.units[index].Ke[2][4] = self.units[index].Ke[4][2] = self.units[index].Ke[2][4] / 2
                        self.units[index].Ke[2][2] = (self.units[index].Ke[2][2] - self.units[index].b) * 3 / 4 + self.units[index].b
                        self.units[index].Ke[1][5] = self.units[index].Ke[2][5] = self.units[index].Ke[4][5] = self.units[index].Ke[5][5] = \
                        self.units[index].Ke[5][4] = self.units[index].Ke[5][2] = self.units[index].Ke[5][1] = 0
                elif len(hitch.points) == 2:
                    for i in [1,2,4,5]:
                        for j in [1,2,4,5]:
                            self.units[index].Ke[i][j] = 0
                return None
            else:
                return None

    def integrateKe(self):
        # calculate global stiffness matrix
        n = len(self.points)
        self.Kg = np.zeros((n*3, n*3))

        for unit in self.units:
            self.calculateHitch(unit)
            unit.transform()
            i, j = unit.i, unit.j
            id = [(i.id-1)*3, i.id*3-2, i.id*3-1, (j.id-1)*3, j.id*3-2, j.id*3-1]

            for p in range(0, 6):
                for q in range(0, 6):
                    self.Kg[id[p]][id[q]] += unit.Ke_ba[p][q]
        # apply constraints
        cl = self.addConstraint()
        mask = [i for i in range(n * 3) if i not in cl]
        self.Kg_calc = self.Kg[np.ix_(mask, mask)]

    def integratePe(self):
        n = len(self.points)
        self.Pe = np.zeros((n*3, 1))
        self.Result = np.zeros(n*3).tolist()
        for p in self.payloads:
            id = p.point.id
            self.Pe[id*3-3][0] += p.px
            self.Pe[id*3-2][0] += p.py
            self.Pe[id*3-1][0] += p.pm
        # apply constraints
        cl = self.addConstraint()
        for i in cl:
            self.Result[i] = "con"
        mask = [i for i in range(n * 3) if i not in cl]
        self.Pe_calc = self.Pe[np.ix_(mask)].T[0]

        return self.Pe_calc
    
    def calculateA(self) -> list[float]:
        inv = np.linalg.inv(self.Kg_calc)
        P = self.Pe_calc.T
        tmp = np.matmul(inv, P).T
        self.A = [np.round(x, 12) for x in tmp]
        for a in self.A:
            if a == 0: # if a is zero, we need to remember the index, in order to replace it again
                a = "con"
            self.Result[self.Result.index(0)] = a
        self.Result = [0 if r == "con" else r for r in self.Result]
        for i, p in enumerate(self.points):
            p.ax = self.Result[i*3]
            p.ay = self.Result[i*3+1]
            p.theta = self.Result[i*3+2] 

        return self.Result
    
    def calculatePe(self):
        self.Pe = np.matmul(self.Kg, self.Result)
        self.Pe = [round(float(x), 3) for x in self.Pe]
        for i, p in enumerate(self.points):
            p.px = self.Pe[i*3]
            p.py = self.Pe[i*3+1]
            p.pm = self.Pe[i*3+2]
        return self.Pe
    
    def transback(self):
        for i, unit in enumerate(self.units):
            T = unit.T
            a = [unit.i.ax, unit.i.ay, unit.i.theta, unit.j.ax, unit.j.ay, unit.j.theta]
            P = np.matmul(unit.Ke, np.matmul(T, a))
            unit.i_N = P[0]
            unit.i_V = P[1]
            unit.i_M = P[2]
            unit.j_N = P[3]
            unit.j_V = P[4]
            unit.j_M = P[5]


class OutputData:
    points:list[Point] = []
    units:list[Unit] = []
    Pe:list[float] = []

    def __init__(self, points:list[Point], units:list[Unit], Pe:list[float]) -> None:
        self.points = points
        self.units = units
        self.Pe = Pe
