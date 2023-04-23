Download Link: https://assignmentchef.com/product/solved-aero4630-project2-conducting-a-tensile-test
<br>
Problem 1: Conducting a Tensile Test

Part 1a: Lam´e Parameters

For <em>E </em>= 200 GPa = 200 × 10<sup>9 </sup>N/m<sup>2 </sup>and <em>ν </em>= 0<em>.</em>3, the Lam´e parameters <em>λ </em>and <em>µ </em>can be calculated as

(1)

and

(2)

Equations 1 and 2 yield <em>λ </em>= 1<em>.</em>153 × 10<sup>11 </sup>N/m<sup>2 </sup>and <em>µ </em>= 7<em>.</em>69 × 10<sup>10 </sup>N/m<sup>2</sup>.

<h2>Part 1b: Tensile Test</h2>

The Paraview output for with stress plotted over the right face of the beam is show in Fig 1

The <em>Slice</em>, <em>ExtractSurface</em>, and <em>IntegrateVariables </em>filters were used to determine that the total force on the face in the x-direction is 936 kN.

<h2>Part 1c: Average Stress vs Strain</h2>

The average stresses and strains for displacements of the right hand side ranging from 0.001 mm to 0.004 mm were plotted, and are shown in Fig 2

The code used for Problem 1 is shown in <strong>Appendix 1</strong>.

Figure 1: Paraview output for problem 1b

Figure 2: Paraview output for problem 1b

<h1>Problem 2: Applying a Force on a Beam Clamped at Both Ends</h1>

<h2>Part 2a: Vertical Displacement given Force Application Area</h2>

A force was applied to an area at the center of beam of length (<em>L</em>−2<em>a</em>) and width <em>w</em>, where <em>a </em>is <em>L/</em>4. The resulting vertical deflection at the middle is −2<em>.</em>322 × 10<sup>−8 </sup>m.

The Python code used for problem 2a is shown in <strong>Appendix 2a</strong>.

<h2>Part 2b: Changing Applied Force Application Area</h2>

Next, displacements at the center of the beam were calculated for different force application areas. The length along the beam at the center at which the force was distributed was <em>a </em>= <em>L/N </em>for <em>N </em>= 3<em>,</em>4<em>,</em>5 and 6. A plot of vertical displacement of the beam over <em>N </em>is shown in Fig 3

Figure 3: Error vs Mesh Sizing for Problem 2a

The code used in part 2b is provided in <strong>Appendix 2b</strong>.

<h1>Appendix 1a: Code for Problem 1</h1>

<em>”””</em>

<em>Python        script      for      Problem 1     of     Project     2</em>

<em>Original          Author :         Vinamra Agrawal</em>

<em>Date :                            January    25 ,    2019</em>

<em>Edited By:                            Omkar Mulekar</em>

<em>Date :                            February     10 ,    2019</em>

<em>”””</em>

<strong>from   </strong>future          <strong> import </strong>print function <strong>from </strong>fenics <strong>import </strong>∗ <strong>from </strong>ufl <strong>import </strong>nabla div

length = .150;

W = 0.025;

H = 0.025;

mu = 7.69 e10

rho = 7800 lambda  = 1.153 e11 g = 10

mesh = BoxMesh( Point (0 ,0 ,0) , Point ( length ,W,H) ,10 ,3 ,3) V = VectorFunctionSpace (mesh , ’P’ ,1) tol = 1E−14

<strong>def </strong>boundary left (x , on boundary) : <strong>return </strong>(onboundary <strong>and </strong>near (x [0] ,0) )

<strong>def </strong>boundary right (x , on boundary) : <strong>return </strong>on boundary <strong>and </strong>near (x [0] , length )

bcleft = DirichletBC (V, Constant ((0 ,0 ,0) ) , boundary left ) bc right = DirichletBC (V, Constant ((0.004 ,0 ,0) ) , boundary right )

<strong>def </strong>epsilon (u) :

<strong>return </strong>0.5∗( nabla grad (u) + nabla grad (u) .T)

<strong>def </strong>sigma(u) :

<strong>return </strong>lambda ∗nabla div (u)∗Identity (d) + mu∗( epsilon (u) + epsilon (u) .T)

u = TrialFunction (V) d = u. geometric dimension () v = TestFunction (V)

f = Constant ((0 ,0 ,0) ) T = Constant ((0 ,0 ,0) )

a = inner (sigma(u) , epsilon (v) )∗dx L = dot( f , v)∗dx + dot (T, v)∗ds

u = Function (V)

solve (a == L, u,                   [ bc left , bc right ])

vtkfileu = File ( ’ deflection . pvd ’ ) vtkfile u <em>&lt;&lt; </em>u

W = TensorFunctionSpace (mesh , ”Lagrange” , 1) stress = lambda ∗nabla div (u)∗Identity (d) + mu∗( epsilon (u) + epsilon (u) .T)

vtkfiles = File ( ’ stress . pvd ’ ) vtkfile s <em>&lt;&lt; </em>project ( stress ,W)

s = sigma(u) − (1./3)∗tr (sigma(u) )∗Identity (d) <em># deviatoric stress </em>von Mises = sqrt (3./2∗ inner (s , s ) )

X = FunctionSpace (mesh , ’P’ , 1) vtkfile von = File ( ’ von Mises . pvd ’ ) vtkfile von <em>&lt;&lt; </em>project ( von Mises , X)

<h1>Appendix 2a: Code for Problem 2a</h1>

<em>”””</em>

<em>Python        script      for      Part 1b     of     Project     2a</em>

<em>Original          Author :         Vinamra Agrawal</em>

<em>Date :                            January    25 ,    2019</em>

<em>Edited By:                            Omkar Mulekar</em>

<em>Date :                            February     10 ,    2019</em>

<em>”””</em>

<strong>from   </strong>future          <strong> import </strong>print function <strong>from </strong>fenics <strong>import </strong>∗ <strong>import </strong>matplotlib

matplotlib . use (”Agg”)

<strong>import </strong>matplotlib . pyplot as plt <strong>from </strong>ufl <strong>import </strong>nabla div

length = .150;

W = 0.025;

H = 0.025;

F = −10 A = length / 4

mu = 7.69 e10

rho = 7800 lambda  = 1.153 e11 g = 10

mesh = BoxMesh( Point (0 ,0 ,0) , Point ( length ,W,H) ,10 ,3 ,3) V = VectorFunctionSpace (mesh , ’P’ ,1) tol = 1E−14

<strong>def </strong>boundary left (x , on boundary) : <strong>return </strong>(onboundary <strong>and </strong>near (x [0] ,0) )

<strong>def </strong>boundary right (x , on boundary) :

<strong>return </strong>on boundary <strong>and </strong>near (x [0] , length )

bcleft = DirichletBC (V, Constant ((0 ,0 ,0) ) , boundaryleft ) bc right = DirichletBC (V, Constant ((0 ,0 ,0) ) , boundary right )

<strong>def </strong>epsilon (u) : <strong>return </strong>0.5∗( nabla grad (u) + nabla grad (u) .T)

<strong>def </strong>sigma(u) :

<strong>return </strong>lambda ∗nabla div (u)∗Identity (d) + mu∗( epsilon (u) + epsilon (u) .T)

u = TrialFunction (V) d = u. geometric dimension () v = TestFunction (V)

f = Constant ((0 ,0 ,0) )

T = Expression (( ’ 0.0 ’ , ’ (x [0] <em>&gt;</em>= a &amp;&amp; x [0] <em>&lt;</em>= (L−a) &amp;&amp; near (x [1] ,W) ) ? (F/(H∗(L−2∗a) ) ) : 0.0 ’ , ’ 0.0 ’ ) ,L=length , a=A,F=F,W=W,H=H, degree=1)

a = inner (sigma(u) , epsilon (v) )∗dx L = dot( f , v)∗dx + dot (T, v)∗ds

u = Function (V)

solve (a == L, u,                   [ bc left , bc right ])

w = u( length /2.0 ,W/2.0 ,H/2.0)

<strong>print</strong>(w[1])

vtkfileu = File ( ’ deflection . pvd ’ ) vtkfile u <em>&lt;&lt; </em>u

W = TensorFunctionSpace (mesh , ”Lagrange” , 1) stress = lambda ∗nabla div (u)∗Identity (d) + mu∗( epsilon (u) + epsilon (u) .T)

vtkfiles = File ( ’ stress . pvd ’ ) vtkfile s <em>&lt;&lt; </em>project ( stress ,W)

s = sigma(u) − (1./3)∗tr (sigma(u) )∗Identity (d) <em># deviatoric stress </em>von Mises = sqrt (3./2∗ inner (s , s ) )

X = FunctionSpace (mesh , ’P’ , 1) vtkfile von = File ( ’ von Mises . pvd ’ ) vtkfile von <em>&lt;&lt; </em>project ( von Mises , X)

<h1>Appendix 2a: Code for Problem 2b</h1>

<em>”””</em>

<em>Python        script      for      Part 1b     of     Project     2b</em>

<em>Original          Author :         Vinamra Agrawal</em>

<em>Date :                            January    25 ,    2019</em>

<em>Edited By:                            Omkar Mulekar</em>

<em>Date :                            February     10 ,    2019</em>

<em>”””</em>

<strong>from   </strong>future          <strong> import </strong>print function <strong>from </strong>fenics <strong>import </strong>∗ <strong>import </strong>matplotlib

matplotlib . use (”Agg”)

<strong>import </strong>matplotlib . pyplot as plt <strong>from </strong>ufl <strong>import </strong>nabla div

length = .150;

W = 0.025;

H = 0.025; F = −10

n = [3 ,4 ,5 ,6]

<em># v e r t i c a l </em><em> d e f l e c t i o n = [ ] </em>∗ <em>len (n) </em>w = [0] ∗ <strong>len</strong>(n)

<h2>for i in range(len(n) ) : print( i )</h2>

A = length / n[ i ]

mu = 7.69 e10

rho = 7800 lambda  = 1.153 e11 g = 10

mesh = BoxMesh( Point (0 ,0 ,0) , Point ( length ,W,H) ,10 ,3 ,3) V = VectorFunctionSpace (mesh , ’P’ ,1) tol = 1E−14

<strong>def </strong>boundary left (x , on boundary) : <strong>return </strong>(onboundary <strong>and </strong>near (x [0] ,0) )

<strong>def </strong>boundary right (x , on boundary) : <strong>return </strong>on boundary <strong>and </strong>near (x [0] , length )

bcleft = DirichletBC (V, Constant ((0 ,0 ,0) ) , boundaryleft ) bc right = DirichletBC (V, Constant ((0 ,0 ,0) ) , boundary right )

<strong>def </strong>epsilon (u) : <strong>return </strong>0.5∗( nabla grad (u) + nabla grad (u) .T)

<strong>def </strong>sigma(u) :

<strong>return </strong>lambda ∗nabla div (u)∗Identity (d) + mu∗( epsilon (u) + epsilon (u) .T)

u = TrialFunction (V) d = u. geometric dimension () v = TestFunction (V)

f = Constant ((0 ,0 ,0) )

T = Expression (( ’ 0.0 ’ , ’ (x [0] <em>&gt;</em>= a &amp;&amp; x [0] <em>&lt;</em>= (L−a) &amp;&amp; near (x

[1] ,W) ) ? (F/(H∗(L−2∗a) ) ) : 0.0 ’ , ’ 0.0 ’ ) ,L=length , a=A,F=F,W=W ,H=H, degree=1)

a = inner (sigma(u) , epsilon (v) )∗dx L = dot( f , v)∗dx + dot (T, v)∗ds

u = Function (V)

solve (a == L, u,            [ bc left , bc right ])

w[ i ] = u( length /2.0 ,W/2.0 ,H/2.0) [1]

<strong>print</strong>(w[1])

<em># v e r t i c a l </em><em> d e f l e c t i o n [ i ] = w[1]</em>

vtkfileu = File ( ’ deflection . pvd ’ ) vtkfile u <em>&lt;&lt; </em>u

W = TensorFunctionSpace (mesh , ”Lagrange” , 1) stress = lambda ∗nabla div (u)∗Identity (d) + mu∗( epsilon (u) + epsilon (u) .T)

vtkfiles = File ( ’ stress . pvd ’ ) vtkfile s <em>&lt;&lt; </em>project ( stress ,W) s = sigma(u) − (1./3)∗tr (sigma(u) )∗Identity (d) <em># deviatoric stress</em>

von Mises = sqrt (3./2∗ inner (s , s ) ) X = FunctionSpace (mesh , ’P’ , 1) vtkfilevon = File ( ’ von Mises . pvd ’ ) vtkfile von <em>&lt;&lt; </em>project ( von Mises , X) W = 0.025

plt . figure (1) plt . plot (n,w, ’b−x ’ ) plt . xlabel ( ’N’ ) plt . ylabel ( ’ Vertical       Deflection       [m] ’ )

plt . savefig ( ’2 afig1 . png ’ )