from . import type, drawer

def writeFile(path: str, output: type.OutputData):

    pic_path = path.replace(".txt", ".png")
    drawer.d(output, pic_path)
    
    try:
        with open(path, "w", encoding="utf8") as file:
            for p in output.points:
                file.write(
                    f"point: {p.id}:({p.x:0.3f}, {p.y:0.3f}),\t displacement: x: {p.ax:0.3e}, y: {p.ay:0.3e}, rotation: {p.theta:0.3e}\n")
            for u in output.units:
                file.write(
                    f"unit: {u.id}, i_N: {u.i_N:0.3e}, i_V: {u.i_V:0.3e}, i_M: {u.i_M:0.3e}, j_N: {u.j_N:0.3e}, j_V: {u.j_V:0.3e}, j_M: {u.j_M:0.3e}\n")
            
            file.write(f"Pe: {output.Pe}")

            # for unit in output.units:
            #     file.write(
            #         f"unit: {unit.id}, strain: {unit.e[0]:0.3e}, {unit.e[1]:0.3e}, {unit.e[2]:0.3e}, stress: {unit.s[0]:0.3e}, {unit.s[1]:0.3e}, {unit.s[2]:0.3e}\n")
    
    except Exception as e:
        print(e, '\n文件输出失败')
