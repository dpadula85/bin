def write_qtiplot_prj(wavelengths, UV, CD, energies, f_osc, rot):
  f = open('prova.qti', 'w')

  f.write('''QtiPlot 0.9.8 project file
<scripting-lang>        muParser
<windows>       2
<table>
Table1  20834   8       16/05/15 11:17
geometry        0       0       445     220
header  lambda[X]    UV spectrum[Y]    CD spectrum[Y]    lambda2[X]    osc[Y]    osc scaled[Y]    rot[Y]    rot scaled[Y]
ColWidth        150     150     150     150     150     150     150     150
<com>
</com>
ColType 0;0/13  0;0/13  0;0/13  0;0/13  0;0/13  0;0/13  0;0/13  0;0/13
ReadOnlyColumn  0       0       0       0       0       0       0       0
HiddenColumn    0       0       0       0       0       0       0       0
Comments
WindowLabel             2
<data>
''')

  for i in range(len(energies)):
    f.write('%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n' % (i, wavelengths[i], UV[i], CD[i], energies[i] * 1240, f_osc[i], f_osc[i] * 50000, rot[i], rot[i] / 3))
    
  for j in range(i + 1, len(wavelengths)):
    f.write('%d\t%f\t%f\t%f\t\t\t\t\t\n' % (j, wavelengths[j], UV[j], CD[j])) 
    
  f.write('''</data>
</table>
<multiLayer>
Graph1  1       1       16/05/15 11:17
geometry        445     0       593     469     active
WindowLabel             2
Margins 5       5       5       5
Spacing 5       5
LayerCanvasSize 400     300
Alignement      0       0
<AlignPolicy>0</AlignPolicy>
<CommonAxes>0</CommonAxes>
<graph>
ggeometry       5       5       575     406
<PageGeometry>0.00854701        0.0118765       0.982906        0.964371</PageGeometry>
PlotTitle       Title   #000000 4228
<Antialiasing>0</Antialiasing>
<Autoscaling>1</Autoscaling>
<ScaleFonts>1</ScaleFonts>
<GridOnTop>0</GridOnTop>
<MissingDataGap>1</MissingDataGap>
Background      #ffffff 0
Margin  0
Border  0       #000000
grid    0       0       0       0       #0000ff 0       0.5     #a0a0a4 2       0.4     #0000ff 0       0.5     #a0a0a4 2       0.4     0       0       2       0       0
EnabledAxes     1       1       1       1
AxesTitles      %(?X)   %(?Y)           %(?Y)
AxesTitleColors #000000 #000000 #000000 #000000
AxesTitleAlignment      5124    5124    5124    5124
AxesTitleDistance       2       2       2       2
InvertedTitle   0       0       0       0
TitleFont       Ubuntu  13      75      0       0       0
ScaleFont0      Ubuntu  11      75      0       0       0
ScaleFont1      Ubuntu  11      75      0       0       0
ScaleFont2      Ubuntu  11      75      0       0       0
ScaleFont3      Ubuntu  11      75      0       0       0
AxisFont0       Ubuntu  11      50      0       0       0
AxisFont1       Ubuntu  11      50      0       0       0
AxisFont2       Ubuntu  11      50      0       0       0
AxisFont3       Ubuntu  11      50      0       0       0
AxesColors      #000000 #000000 #000000 #000000
AxesNumberColors        #3c3c3c #3c3c3c #3c3c3c #3c3c3c
AxesBaseline    0       0       0       0
CanvasBackground        #ffffff 0
curve   Table1_1        Table1_3        2       1       #000000 0       1       7       1       #000000 #000000 0       #000000 14      1       2       0       0       20833   1
curve   Table1_1        Table1_2        2       1       #ff0000 0       1       7       2       #ff0000 #ff0000 0       #000000 14      1       2       1       0       20833   1
scale   0       -600    600     0       7       5       0       0
scale   1       -50000  350000  0       9       5       0       0
scale   2       100     450     0       8       5       0       0
scale   3       100     450     0       8       5       0       0
LabelsFormat    0       6       0       6       0       6       0       6
AxisType        0       0       0       0
MajorTicks      1       1       1       1
MinorTicks      1       1       1       1
TicksLength     5       9
DrawAxesBackbone        1       1       1       1       1
AxesLineWidth   1
LabelsRotation  0       0       0       0
LabelsPrefix
LabelsSuffix
TickLabelsSpace 4       4       4       4
ShowTicksPolicy 0       0       0       0
EnabledTickLabels       1       1       1       1
<Legend>
<Frame>1</Frame>
<Color>#000000</Color>
<FrameWidth>1</FrameWidth>
<LineStyle>0</LineStyle>
<x>108.77192982456</x>
<y>336622.0735786</y>
<right>195.614035087719</right>
<bottom>297826.086956522</bottom>
<attachTo>1</attachTo>
<onTop>1</onTop>
<visible>1</visible>
<Text>
\l(1)%(1)
\l(2)%(2)
</Text>
<Font>Ubuntu    11      50      0       0       0</Font>
<TextColor>#000000</TextColor>
<Background>#ffffff</Background>
<Alpha>0</Alpha>
<Angle>0</Angle>
<AutoUpdate>1</AutoUpdate>
<TeXOutput>0</TeXOutput>
</Legend>
</graph>
<LinkXAxes>0</LinkXAxes>
<ScaleLayers>1</ScaleLayers>
</multiLayer>
<open>1</open>
''')