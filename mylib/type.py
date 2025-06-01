import numpy as np
import math

class Point: # 节点
    id = 0
    # coordinate
    x = 0
    y = 0
    theta = 0
    # load
    px = 0
    py = 0
    pm = 0
    # displacement
    ax = 0
    ay = 0
    # rotation
    a_theta = 0

    def __init__(self, data:list[float]):
        (self.id, self.x, self.y) = data

    def __str__(self):
        return f"Node {self.id}: ({self.x}, {self.y})\n"
    
    def __repr__(self):
        return self.__str__()
    
    def update(self) -> None:
        # update coordinate
        self.x += self.ax
        self.y += self.ay
        self.theta += self.a_theta
        # self.theta += self.a_theta
        self.ax = 0
        self.ay = 0
        self.a_theta = 0
        return None

class Material: # 材料
    id = 0
    E = 0
    v = 0
    G = 0

    def __init__(self, data: list[float]) -> None:
        (self.id, self.E, self.v) = data

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
        # calculate moment of inertia
        if self.shape == 0: # circle
            self.A = math.pi * self.b**2 / 4
            self.I = math.pi * self.b**4 / 64
        elif self.shape == 1: # rectangle
            self.A = self.b * self.h
            self.I = self.b * self.h**3 / 12

    def __str__(self) -> str:
        return f"<id:{self.id}, A:{self.A}, I:{self.I}, shape:{self.shape}>"

    def __repr__(self) -> str:
        return self.__str__()

class Unit: # 单元
    #region parameters
    id: int = 0
    # node id
    i: Point = None
    j: Point = None
    # material
    material: Material = None
    # plane
    plane: Plane = None
    # length
    L = 0
    # rigid angle
    theta = 0
    # trsnformation matrix
    T = np.zeros((6, 6))
    # axial stiffness
    EA = 0
    # flexural stiffness
    EI = 0
    # force vector
    unit_F:np.array = np.zeros((6,))  # [N1, V1, M1, N2, V2, M2]
    unit_u:np.array = np.zeros((6,))  # [u1, v1, theta1, u2, v2, theta2]
    delta_u:np.array = np.zeros((6,))  # [0, 0, i_angle, u, 0, j_angle]
    # endregion

    def __init__(self, data:list[float]) -> None:
        (self.id, self.i, self.j, self.material, self.plane) = data
        # calculate length
        dx = self.j.x - self.i.x
        dy = self.j.y - self.i.y
        self.L = np.hypot(dx, dy)
        # calculate angle
        self.theta = np.arctan2(dy, dx)
        # calculate axial stiffness
        self.EA = self.material.E * self.plane.A
        # calculate flexural stiffness
        self.EI = self.material.E * self.plane.I
        self.unit_F = np.zeros((6,))  # [N1, V1, M1, N2, V2, M2]
        self.unit_u = np.zeros((6,))  # [u1, v1, theta1, u2, v2, theta2]
        self.delta_u = None  # [0, 0, i_angle, u, 0, j_angle]

    def __str__(self):
        return f"Element {self.id}: ({self.i}, {self.j})"

    def __repr__(self):
        return self.__str__()
    
    def update(self):
        # update length
        dx = self.j.x - self.i.x
        dy = self.j.y - self.i.y
        self.L = np.hypot(dx, dy)
        # update angle
        self.theta = np.arctan2(dy, dx)
        # update transformation matrix
        self.calculateT()
    
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
        # self.update()  # update length and angle
        self.calculateT()  # update transformation matrix
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
        self.Ke = Ke
        return Ke
    
    def calculateStressKe(self) -> np.ndarray:
        # element_info:总信息中的单元信息
        # alpha:顺时针为正
        # element_force: 单元坐标系下单元内力，[N1, V1, M1, N2, V2, M2]
        # self.update()  # update length and angle
        I = self.plane.I
        A = self.plane.A
        l = self.L
        
        Fx = self.unit_F[3]
        Mi = self.unit_F[2]
        Mj = self.unit_F[5]
        FI = Fx * I
        FL = Fx * l
        AL = A * l
        AL2 = A * l**2
        AL3 = A * l**3
        Kg = np.array([[Fx/l,          0,                           -Mi/l,           -Fx/l,         0,                          -Mj/l],
                        [0,         12*FI/AL3 + 6*Fx/(5*l),      6*FI/AL2 + Fx/10,      0,        -12*FI/AL3 - 6*Fx/(5*l),    6*FI/AL2 + Fx/10],
                        [-Mi/l,     6*FI/AL2 + Fx/10,            4*FI/AL + 2*FL/15,    Mi/l,      -6*FI/AL2 - Fx/10,          2*FI/AL - FL/30],
                        [-Fx/l,         0,                            Mi/l,            Fx/l,          0,                           Mj/l],
                        [0,         -12*FI/AL3 - 6*Fx/(5*l),     -6*FI/AL2 - Fx/10,     0,        12*FI/AL3 + 6*Fx/(5*l),     -6*FI/AL2 - Fx/10],
                        [-Mj/l,     6*FI/AL2 + Fx/10,            2*FI/AL - FL/30,      Mj/l,      -6*FI/AL2 - Fx/10,          4*FI/AL + 2*FL/15]])
        # Kg = np.dot(self.T.T, Kg)
        # Kg = np.dot(Kg, self.T)
        self.Kg = Kg
        return Kg

    def split_displacement(self):
        self.unit_u[0] = self.i.ax
        self.unit_u[1] = self.i.ay
        self.unit_u[2] = self.i.a_theta
        self.unit_u[3] = self.j.ax
        self.unit_u[4] = self.j.ay
        self.unit_u[5] = self.j.a_theta
        self.unit_u = np.dot(self.T, self.unit_u)
        l = self.L
        dx = l + self.unit_u[3] - self.unit_u[0]
        dy = self.unit_u[4] - self.unit_u[1]
        l2 = np.hypot(dx, dy)
        u = l2 - l
        rigid_angle = np.arctan2(dy, dx)
        i_angle = self.unit_u[2] - rigid_angle
        j_angle = self.unit_u[5] - rigid_angle
        delta_u = np.array([0, 0, i_angle, u, 0, j_angle])
        self.delta_u = delta_u
        return delta_u
    
    def calculateF(self):
        K = self.calculateKe_without_shear() + self.calculateStressKe()
        u = self.delta_u
        F = np.matmul(K, u)
        return F


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

class GeoNonlinear: # 几何非线性
    points: list[Point] = []
    planes: list[Plane] = []
    units: list[Unit] = []
    materials: list[Material] = []
    payloads: list[Payload] = []
    constraints: list[Constraint] = []

    def __init__(self) -> None:
        self.points, self.planes, self.units, self.materials, self.payloads, self.constraints = [], [], [], [], [], []

    def __str__(self) -> str:
        return f"<\npoints:{self.points},\nplanes:{self.planes},\nunits:{self.units},\nmaterials:{self.materials},\npayloads:{self.payloads},\nconstraints:{self.constraints}\n>"

    def __repr__(self) -> str:
        return self.__str__()

    def reset(self) -> None:
        self.Result = np.zeros(self.points_num*3).tolist()
        cl = self.cl
        for i in cl:
            self.Result[i] = "con"

    def calculateNumber(self) -> None:
        self.points_num = len(self.points)
        self.units_num = len(self.units)
        return None

    def addConstraint(self) -> list[int]:
        ConstraintLine = []
        for constraint in self.constraints:
            ConstraintLine.append(int((constraint.point.id-1)*3+constraint.axis-1))
        self.cl = ConstraintLine
        return ConstraintLine

    def integrateKe(self) -> np.ndarray:
        Kg = np.zeros((self.points_num * 3, self.points_num * 3))    # 全局刚度矩阵
        for unit in self.units:                                     # 注意更新节点坐标，此事为t时刻坐标
            K_L = unit.calculateKe_without_shear()
            K_NL = unit.calculateStressKe()
            K = K_L + K_NL
            K = np.dot(np.dot(unit.T.T, K), unit.T)
            i, j = unit.i, unit.j
            id = [(i.id-1)*3, i.id*3-2, i.id*3-1, (j.id-1)*3, j.id*3-2, j.id*3-1]

            for p in range(0, 6):
                for q in range(0, 6):
                    Kg[id[p]][id[q]] += K[p][q]

        # apply constraints
        cl = self.addConstraint()
        mask = [i for i in range(self.points_num * 3) if i not in cl]
        Kg_calc = Kg[np.ix_(mask, mask)]

        self.Kg = Kg
        self.Kg_calc = Kg_calc
        return Kg_calc
    
    def integratePe(self):
        Pe = np.zeros((self.points_num*3, 1))
        Result = np.zeros(self.points_num*3).tolist()
        for p in self.payloads:
            id = p.point.id
            Pe[id*3-3][0] += p.px
            Pe[id*3-2][0] += p.py
            Pe[id*3-1][0] += p.pm

        # apply constraints
        cl = self.addConstraint()
        for i in cl:
            Result[i] = "con"
        mask = [i for i in range(self.points_num * 3) if i not in cl]
        Pe_calc = Pe[np.ix_(mask)].T[0]

        self.Pe = Pe
        self.Pe_calc = Pe_calc
        self.Result = Result

        return Pe_calc
    
    def calculateA(self, P: np.array) -> list[float]:
        inv = np.linalg.inv(self.Kg_calc)
        # P = self.Pe_calc.T
        cl = self.addConstraint()
        mask = [i for i in range(self.points_num * 3) if i not in cl]
        P = P[np.ix_(mask)].T[0]

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
            p.a_theta = self.Result[i*3+2] 

        return self.Result

    def integrateF(self):
        internal_F = np.zeros((self.points_num * 3, 1))
        for unit in self.units:
            tmp = np.matmul(unit.T.T, unit.unit_F)
            i, j = unit.i, unit.j
            id = [(i.id-1)*3, i.id*3-2, i.id*3-1, (j.id-1)*3, j.id*3-2, j.id*3-1]
            for p in range(0, 6):
                internal_F[id[p]][0] += tmp[p]
        self.internal_F = internal_F
        return internal_F
        
    def iterate(self, steps: int=200, max_iterator:int = 10, error: float=1e-3) -> None:
        # 迭代求解
        self.calculateNumber()
        self.integratePe()
        F = [np.zeros((self.points_num * 3, 1))]  # 初始化F
        u = [np.zeros((self.points_num * 3, 1))]
        outputs:list[OutputData] = []
        
        for step in range(steps):
            P_step = self.Pe * ((1 + step) / steps)  # 逐步加载
            F_step = F[step]
            u_step = u[step]
            R_step = P_step - F_step
            for iter in range(max_iterator):    # 迭代求解
                self.reset()
                if max(abs(R_step)) >= error * max(abs(P_step)):
                    print(f"step: {step}, iteration: {iter}")
                    self.integrateKe()
                    self.calculateA(R_step)
                    for unit in self.units:
                        unit.split_displacement()
                        delta_F = unit.calculateF()
                        unit.unit_F += delta_F
                    F_step = self.integrateF()
                    u_step += np.array(self.Result).reshape((self.points_num * 3, 1))

                    for point in self.points:
                        point.update()
                    for unit in self.units:
                        unit.update()
                    R_step = P_step - F_step
                    for c in self.cl:
                        R_step[c] = 0
                        
                    if iter == max_iterator - 1:
                        print(f"Not converged in step {step}, iteration {iter}")
                        break
                else:
                    outputs.append(OutputData(self.points).__copy__())
                    F.append(F_step)
                    u.append(u_step)
                    break

        return outputs
            
class OutputData:
    points:list[Point] = []

    def __init__(self, points:list[Point]) -> None:
        self.points = points

    def __copy__(self):
        new_points = [Point([p.id, p.x, p.y]) for p in self.points]
        return OutputData(new_points)