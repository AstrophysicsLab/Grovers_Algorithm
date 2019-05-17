N = 7;
Question = 7;
Measurments =1E5;
Psi = Qbit();
Psi.N = N;
Psi.Initialize_Zero_State(); 
Psi.Generate_Grover_Block(Question)
Psi.Grovers_Algorithm();
Psi.Measure(Measurments);
%Psi.Measurement_Array
%  Psi.State_Matrix
Psi.Measured_Ket
%  Psi.Measured_Matrix
%  Psi.Bar_Chart_Results();
% Psi.Normalized_Measurement