cd %~dp0

triangle.exe -pq33aDenV boundtriangle
triangle2plt.exe boundtriangle.1 su2
triangle2plt.exe boundtriangle.1 plt


mv boundtriangle.1.su2 triangularmesh.su2
mv boundtriangle.1.plt triangularmesh.plt

find . -maxdepth 1 -not -name "triangularmesh*" -delete 
