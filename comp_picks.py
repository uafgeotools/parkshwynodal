import os
def compare_directories(dir1, dir2, dir3):
    x = 0
    for root, _, files in os.walk(dir1):
        for filename in files:
            print(root)
            file1 = os.path.join(root, filename)
            file2 = os.path.join(dir2, os.path.relpath(file1, dir1))
            file3 = os.path.join(dir3, os.path.relpath(file1, dir1))
            x += 3
            if not os.path.exists(file2):
                print(f"Missing in {dir2}: {file2}")
            if not os.path.exists(file3):
                print(f"Missing in {dir3}: {file3}")
    return x
yy = 0
for dirc in os.listdir('input/Data_Picks'):
    dirc_list = []
    new_dirc = os.path.join('input/Data_Picks', dirc)
    for dirc in os.listdir(new_dirc):
        dirc_list.append(os.path.join(new_dirc, dirc))
    xx = compare_directories(dirc_list[0], dirc_list[1], dirc_list[2])
    yy += xx
print(yy)
