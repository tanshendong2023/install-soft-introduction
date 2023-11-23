### 一、建立模型


```python
# 输入参数
num_molecules1 = 20  # 第一种分子的数量
num_molecules2 = 20  # 第二种分子的数量
density = 3.0        # mol/L -r
distance = 2.5       # 分子间距 -t
file_path1 = "./data/c4c1im.zmat"  # 第一个文件的路径
file_path2 = "./data/ntf2.zmat"   # 第二个文件的路径

# 构建命令 build a simulation box with x and y molecules and a density of densiy mol/L
command = f"fftool {num_molecules1} {file_path1} {num_molecules2} {file_path2} -r {density}"

# 在 Jupyter Notebook 中执行命令fftool
!{command}

# 用packmol建模
import os

# Add the path of Packmol
#os.environ['PATH'] -= os.pathsep + '/opt/gengzi/pub/softwares/Packmol/packmol-20.3.5'

# Use packmol with the pack.inp file
!packmol < pack.inp >> output.txt
```

    density 3.000 mol/L  volume 22141.0 A^3
    molecule_file      species           nmol force_field      nbond source  charge
      ./data/c4c1im.zmat c4c1im+             20 il.ff               25 file   +1.0000
      ./data/ntf2.zmat tf2N-               20 il.ff               14 file   -1.0000
    packmol file
      pack.inp


### 二、可视化模型


```python
import subprocess
import os
import tempfile

def format_and_convert_zmatrix(input_file_path, xyz_file_path):
    # 读取并格式化 Z-matrix
    with open(input_file_path, 'r') as file:
        lines = file.readlines()

    formatted_lines = []
    for line in lines:
        parts = line.split()
        if parts and parts[0].isdigit():
            formatted_line = parts[1][0] + ' ' + ' '.join(parts[2:])
            formatted_lines.append(formatted_line)

    # 将格式化后的数据保存到临时文件中
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.zmat') as temp_file:
        temp_file.writelines('\n'.join(formatted_lines))
        temp_file_path = temp_file.name

    # 转换格式化后的 Z-matrix 到 XYZ 格式
    with open(xyz_file_path, 'w') as xyz_file:
        subprocess.run(
            ["python", "/home/tanshendong/soft/geomConvert-master/gc.py", "-zmat", temp_file_path],
            stdout=xyz_file
        )

    # 删除临时文件
    os.remove(temp_file_path)

# 输入文件和 XYZ 文件路径
input_file_path_1 = './data/c4c1im.zmat'
input_file_path_2 = './data/ntf2.zmat'
xyz_file_path_1 = 'new_c4c1im.xyz'
xyz_file_path_2 = 'new_ntf2.xyz'

# 调用函数处理文件并转换格式
format_and_convert_zmatrix(input_file_path_1, xyz_file_path_1)
format_and_convert_zmatrix(input_file_path_2, xyz_file_path_2)
print("Conversion to XYZ format is complete.")
```

    Conversion to XYZ format is complete.



```python
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from IPython.display import display
import ipywidgets as widgets
import io

def visualize_molecule_from_xyz(xyz_filename):
    # 使用 Open Babel 从 XYZ 文件读取分子
    mol = next(pybel.readfile("xyz", xyz_filename))

    # 将分子转换为 SMILES
    smiles = mol.write("smi")

    # 使用 RDKit 从 SMILES 创建分子
    rdkit_mol = Chem.MolFromSmiles(smiles.strip())

    # 对分子进行二维坐标生成
    rdkit_mol = Chem.AddHs(rdkit_mol)
    AllChem.Compute2DCoords(rdkit_mol)

    # 可视化分子
    return Draw.MolToImage(rdkit_mol)

# 可视化两个分子
img1 = visualize_molecule_from_xyz('new_c4c1im.xyz')
img2 = visualize_molecule_from_xyz('new_ntf2.xyz')

# 将图像转换为可显示的格式
bio1 = io.BytesIO()
bio2 = io.BytesIO()
img1.save(bio1, format='PNG')
img2.save(bio2, format='PNG')

image_widget1 = widgets.Image(value=bio1.getvalue(), format='png', width=200, height=200)
image_widget2 = widgets.Image(value=bio2.getvalue(), format='png', width=200, height=200)

# 使用 HBox 布局并排显示
hbox = widgets.HBox([image_widget1, image_widget2])
display(hbox)
```

    ==============================
    *** Open Babel Warning  in PerceiveBondOrders
      Failed to kekulize aromatic bonds in OBMol::PerceiveBondOrders (title is INSERT TITLE CARD HERE)
    



    HBox(children=(Image(value=b'\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x01,\x00\x00\x01,\x08\x02\x00\x00\x00…



```python
import py3Dmol

# 读取.xyz文件
with open('./simbox.xyz', 'r') as file:
    xyz_data = file.read()

# 设置视图
# 创建视图
view = py3Dmol.view(width=800, height=400)
# 添加分子模型
view.addModel(xyz_data, "xyz")
# 设置分子的显示风格
view.setStyle({'stick': {}})
# 调整视图
view.zoomTo()
# 显示视图
view.show()


### 三、生成GROMACS输入文件


```python
# Use fftool to create input files for GROMACS (-g)
command = f"fftool {num_molecules1} {file_path1} {num_molecules2} {file_path2} -r {density} -g"

# 在 Jupyter Notebook 中执行命令fftool
!{command}
```

    density 3.000 mol/L  volume 22141.0 A^3
    molecule_file      species           nmol force_field      nbond source  charge
      ./data/c4c1im.zmat c4c1im+             20 il.ff               25 file   +1.0000
      ./data/ntf2.zmat tf2N-               20 il.ff               14 file   -1.0000
    gromacs files
      em.mdp
      npt-gas.mdp
      npt.mdp
      field.top
      config.pdb


### 四、运行MD模拟基于GROMACS


```python
import subprocess

# 定义一个函数来执行命令并捕获输出
def run_command(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    return stdout.decode(), stderr.decode()

# 执行 gmx_mpi editconf
stdout, stderr = run_command(['gmx', 'editconf', '-f', 'config.pdb', '-o', 'conf.gro'])

# 能量极小化
# 执行 gmx_mpi grompp
stdout, stderr = run_command(['gmx', 'grompp', '-f', 'em.mdp', '-c', 'conf.gro', '-p', 'field.top', '-o', 'em'])
# 执行 gmx_mpi mdrun
stdout, stderr = run_command(['gmx', 'mdrun', '-v', '-deffnm', 'em'])

# npt模拟 10ns
# 执行 gmx_mpi grompp
stdout, stderr = run_command(['gmx', 'grompp', '-f', 'npt.mdp', '-c', 'em.gro', '-p', 'field.top', '-o', 'npt'])
# 执行 gmx_mpi mdrun, 提交任务到超算中心
!sbatch sub.gromacs-cpu
```

    Submitted batch job 112635


### 五、MD模拟后处理分析


```python
### 待续
```
