from . import type, drawer
import re


def readLine(line: str) -> str:
    tmp = line.split('#')
    return re.sub('[\\s,，；;]+', ' ', tmp[0]).strip()

def readFile(filePath: str, shear: bool = False) -> type.InputData:
    fileData = []
    pic_path = filePath.replace(".txt", ".png")
    try:
        with open(filePath, mode='r', encoding='utf8') as file:
            for i in file:
                tmp = readLine(i)
                if (tmp == ''):
                    continue
                tmp = [float(x)
                       for x in filter(lambda x: x != '', tmp.split())]
                tmp[0] = int(tmp[0])
                fileData.append(tmp)
    except Exception as e:
        print(e, '\n文件读取失败，请检查文件路径')
    return parseData(fileData, pic_path, shear)

def parseData(src: list[list[float]], pic_path: str, shear: bool = False) -> type.InputData:
    data = type.InputData()
    try:
        # 节点数目和节点坐标：编号  x坐标  y坐标
        count = src.pop(0)[0]
        while count > 0:
            data.points.append(type.Point(src.pop(0)))
            count -= 1
        # 材料信息 ：弹性模量  泊松比
        count = src.pop(0)[0]
        while count > 0:
            data.materials.append(type.Material(src.pop(0)))
            count -= 1
        # 截面信息：编号  厚度
        count = src.pop(0)[0]
        while count > 0:
            data.planes.append(type.Plane(src.pop(0)))
            count -= 1
        # 单元信息：单元编号 i点 j点 m点 材料号  截面号
        count = src.pop(0)[0]
        while count > 0:
            tmp = src.pop(0)
            tmp[1] = getItemById(data.points, tmp[1])
            tmp[2] = getItemById(data.points, tmp[2])
            tmp[3] = getItemById(data.materials, tmp[3])
            tmp[4] = getItemById(data.planes, tmp[4])
            data.units.append(type.Unit(tmp, shear))
            count -= 1
        # 荷载信息 ：节点号，x方向力，y方向力
        count = src.pop(0)[0]
        while count > 0:
            tmp = src.pop(0)
            tmp[0] = getItemById(data.points, tmp[0])
            data.payloads.append(type.Payload(tmp))
            count -= 1
        # 位移约束信息：节点号，约束方向，约束值
        count = src.pop(0)[0]
        while count > 0:
            tmp = src.pop(0)
            tmp[0] = getItemById(data.points, tmp[0])
            data.constraints.append(type.Constraint(tmp))
            count -= 1
        # 铰接: 单元号，节点号
        count = src.pop(0)[0]
        while count > 0:
            tmp = src.pop(0)
            tmp[0] = getItemById(data.units, tmp[0])
            tmp[1] = getItemById(data.points, tmp[1])
            if data.hitches == []:
                data.hitches.append(type.Hitch(tmp[0]))
                data.hitches[-1].addPoint(tmp[1])
            else:
                for hitch in data.hitches:
                    if hitch.unit == tmp[0]:
                        index = data.hitches.index(hitch)
                        data.hitches[index].addPoint(tmp[1])
                    else:
                        data.hitches.append(type.Hitch(tmp[0]))
                        data.hitches[-1].addPoint(tmp[1])
            count -= 1

    except Exception as e:
        Exception(e, '\n输入的文件格式不正确，数据文件出错')

    drawer.d(data, pic_path)
    return data

def getItemById(list, id):
    for i in list:
        if i.id == id:
            return i
    Exception(f"{list}中没有指定的资源编号：{id}，数据文件出错")

if __name__ == "__main__":
    data = readFile("..\\input\\1.txt")
    print(data)