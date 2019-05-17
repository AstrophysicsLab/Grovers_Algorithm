classdef Qbit < handle
    %Qbit 
    %   Class to represent Quantum State
    
    properties
        A;
        N;
        n;
        N_Binary_Ket;
        State_Ket;
        Initial_Ket;
        Initial_State;
        State_Matrix;
        Measured_Ket;
        Measured_Matrix;
        Measurement_Array = [];
        Normalized_Measurement;
        H = [1,1;1,-1];
        H_n;
        H_Array;
        H_n_Sparse;
        Hadamard_Values;
        R1;
        R2;
        R3;
        O;
        J;
        q;
        q_prime;
        Qbits;
        Array_index = [0,0,0];
        Indices;
        Random_Array;
    end
    
    methods

       
       %%%%%%%%%%%%%%%% Initilize Class %%%%%%%%%%%%%%%%%%%%%%%
        
       function Rand_Basis_Ket = Random_Basis_Ket(obj)
           Rand_Basis_Ket = randi([0 1],1,3);    
       end
       
       function Initialize_Zero_State(obj) 
           obj.State_Matrix = obj.Qbit_Basis_State(1);
           obj.Measurement_Array = zeros(2^obj.N,1);
           obj.Initial_State = obj.State_Matrix;
           obj.State_Ket = Build_N_Qbit_Vector(obj,1);
       end
       
       
        function  Matrix = Qbit_Basis_State(obj, idx)
           Matrix = zeros(2^obj.N,1);
           Matrix(idx,1) = 1;
        end
       
        
        function Random_State = Initialize_Random_Basis_State(obj)
             r = randi([1 2^obj.N],1,1);
             Random_State = zeros(2^obj.N,1);
             Random_State(r)= 1;
         end
       
        %%%%%%%%%%%%%%% Class Utilities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        function Ket= Build_State_Matrix(obj,idx)
            index = find(Matrix);
            Ket = obj.Build_N_Qbit_Vector(index);
        end
        
        function Normalize_State(obj)
            obj.Normalized_Measurement =(1/sqrt(sum(obj.Measured_Matrix)))*obj.State_Matrix;
        end
        
          function Binary_Vector = Build_N_Qbit_Vector(obj,index)
            Binary_Vector = zeros(1,obj.N);
            Index_Binary_Decimal = dec2bin(index-1);
            start = obj.N - length(Index_Binary_Decimal);
            for idx =  1:length(Index_Binary_Decimal)
              Binary_Vector(start+idx) = str2double(Index_Binary_Decimal(idx));
            end
         end
         
         function Generate_Qbit_Binary_idx(obj)
             obj.Qbits = zeros(2^obj.N,obj.N);
             for idx = 1:(2^obj.N)
                 obj.Qbits(idx,:) = obj.Build_N_Qbit_Vector(idx);
             end
         end
       
         function Random_State = Random_Basis_State(obj, idx)
             Random_State = zeros(2^obj.N,1);
             Random_State(obj.Random_Array(idx))= 1;
         end
       
       %%%%%%%%%%%%%%%% Measurments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       

       
       function Probability = Probability(obj,Bra)
            Probability = abs(Bra'*obj.State_Matrix)^2;
       end
       
       function [Probability, index] =Random_Measurement(obj,idx)     
         Random_State = obj.Random_Basis_State(idx);
         index = find(Random_State);
         Probability = obj.Probability(Random_State);       
       end
       
       function Parse_Random_Measurement(obj,Probability, Matrix_Value)
           if Probability > rand
                obj.Measurement_Array(Matrix_Value) = obj.Measurement_Array(Matrix_Value)+1;
           end
       end
       
  
       
       function Parse_Results(obj)
           [~,idx] = max(obj.Measurement_Array);
           obj.Measured_Ket = obj.Build_N_Qbit_Vector(idx); 
           obj.Measured_Matrix = obj.Qbit_Basis_State(idx);
       end
       
       function Measure(obj, Measurements)
                       
           obj.Random_Array = randi([1 2^obj.N],Measurements,1);

           for idx = 1:Measurements
              [Probability, Matrix_Value] =obj.Random_Measurement(idx); 
              obj.Parse_Random_Measurement(Probability, Matrix_Value);
           end
%             while sum(obj.Measurement_Array)< N
%                  [Probability, Matrix_Value] =obj.Random_Measurement(); 
%                   obj.Parse_Random_Measurement(Probability, Matrix_Value);
%             end
           obj.Parse_Results();
           obj.Normalize_State();
       end
           
 
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%% Matrices and Gates%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       %%%%%%%%%%%%% Oracle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       function Generate_DxD_Oracle(obj,Question)
              obj.O = eye(2^obj.N,2^obj.N);
              obj.O(Question,Question) = -1;
           end


       %%%%%%%%%%%%%Hadarmard%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       function Generate_All_Hadamard_Gates(obj)
               obj.Generate_Qbit_Binary_idx();
               obj.H_Array = zeros(2^obj.N,2^obj.N,obj.N);
               for idx = 1:obj.N
                   obj.Kronecker_delta(idx);
                   obj.Find_Q(idx);
                   obj.Get_Hadamard_Values();
                    for jdx = 1:length(obj.Indices)
                         obj.H_Array(obj.Indices(jdx,1),obj.Indices(jdx,2),idx) =obj.Hadamard_Values(jdx);
                    end
                    
               end
               
        end
           
        function Apply_Hadamard(obj,Hn)
           obj.State_Matrix =Hn*obj.State_Matrix;
       end
       
 
       
        function Get_Hadamard_Values(obj)
                obj.Hadamard_Values = zeros(length(obj.Indices),1);
                for idx = 1:length(obj.Indices)
                 obj.Hadamard_Values(idx) =obj.H(obj.q(obj.Indices(idx,1)), obj.q(obj.Indices(idx,2)));
                end
                obj.Hadamard_Values =  1/sqrt(2)*obj.Hadamard_Values;
                
        end
        
        
        function Generate_Hadamard_Matrix(obj)
                 obj.H_n = zeros(2^obj.N,2^obj.N);
                 for idx = 1:length(obj.Indices)
                    obj.H_n(obj.Indices(idx,1),obj.Indices(idx,2)) =obj.Hadamard_Values(idx);
                 end
                 
        end
        
        
        function Generate_Sparse_Hadamard_Matrix(obj)
   
              obj.H_n_Sparse = sparse(obj.Indices(:,1),obj.Indices(:,2),obj.Hadamard_Values,2^obj.N,2^obj.N); 
        end
        
        function Generate_Hadamard_Gate(obj,n)
               obj.Generate_Qbit_Binary_idx();
               obj.Kronecker_delta(n);
               obj.Find_Q(n);
               obj.Get_Hadamard_Values();
               obj.Generate_Hadamard_Matrix();
         
        end
        
        %%%%%%%%%%%%%%% Matrice Utilities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               
%         function Kronecker_delta(obj,n)
%              obj.Array_index =zeros(2^obj.N*2,3);
%              Temp_Qbit_index = obj.Qbits;
%              Temp_Qbit_index(:,n) = [];
%              position = 1;
%              index = [0,0];
%            
%              for idx = 1:2^obj.N
%                  index(1,1) = idx;
%                  for jdx = 1: 2^obj.N
%                      index(1,2) = jdx;
%                       if isequal(Temp_Qbit_index(idx,:),Temp_Qbit_index(jdx,:))
%                         obj.Indices(position,:) = index;
%                         position = position+1;
%                       end                 
%                  end
%              end    
%         end
%         
            

  
         function Kronecker_delta(obj,n)
            obj.Array_index =zeros(2^obj.N*2,3);
            Temp_Qbit_index = obj.Qbits;
            Temp_Qbit_index(:,n) = [];
            String_Qbit_index =string(Temp_Qbit_index);
            [row,col] = size(String_Qbit_index);
            Decimal_Qbit_index = zeros(row,1);
            for idx = 1:col 
                Decimal_Qbit_index(:) = Decimal_Qbit_index(:)+ String_Qbit_index(:,idx);
            end
            
            obj.A = Decimal_Qbit_index;
            %Decimal_Qbit_index = str2double(Decimal_Qbit_index); 
            position = 1;
            index = [0,0]; 
            for idx = 1:2^obj.N 
                found = 0;
                index(1,1) = idx;
                for jdx = idx: 2^obj.N 
                    index(1,2) = jdx; 
                    
                        if isequal(Decimal_Qbit_index(idx),Decimal_Qbit_index(jdx))  
                            obj.Indices(position,:) = index;
                            position = position+1;
                            found = found + 1;
                        end
                    
                 end
            end    
             
             mat = zeros(2^obj.N,2^obj.N);
             len = 2*2^obj.N- (1/2*(2^obj.N));
             for idx =1:len
                 mat(obj.Indices(idx,1),obj.Indices(idx,2)) = 1;
             end
   
             mat = mat'+mat;
             obj.Indices = zeros(2^obj.N*2,2);
             
             [obj.Indices(:,1),obj.Indices(:,2)] = find(mat); 
        end
        
        function Find_Q(obj,n)
            obj.q = zeros(2^obj.N,1);
            for idx = 1:2^obj.N
                obj.q(idx) = obj.Qbits(idx,n)+1;
            end
         end

       %%%%%%%%%%%%%%%%Phase Shift%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function Generate_Phase_Shift_Gates(obj)
          %Generate_Phase_Shift_Gates
          % Generates the Phase Shift Gates
          Euler = exp(1i*pi);
          R = [1,0 ; 0, Euler];
          I = [1, 0; 0, -1];
          obj.R1 =   kron(R,kron(I,I));
          obj.R2 =  kron(I,kron(R,I));
          obj.R3 =  kron(I,kron(I,R));
       end

       function Generate_J(obj)
           obj.J = eye(2^obj.N,2^obj.N);
           obj.J(1,1) = -1;    
       end
             
          %%%%%%%%%%%%%%%%%%%Apply Gates%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
          function Apply_All_Hadamard(obj)
              % Apply_All_Hadamard
              % Applies the Hadamard Gate of every Qbit to the state
             for idx = 1:obj.N
                obj.State_Matrix =obj.H_Array(:,:,idx)*obj.State_Matrix;
             end
           end
        
        function Apply_Phase_Shift(obj,Rn)
           obj.State_Matrix =Rn*obj.State_Matrix;
        end
       
        function Apply_Oracle(obj)
            % Apply_Oracle
            % Applies the Oracle to the state
            obj.State_Matrix =obj.O*obj.State_Matrix;
        end
        
         function Apply_J(obj)
            obj.State_Matrix =obj.J*obj.State_Matrix;
         end
         
         
         %%%%%%%%%%%%%%%%%%%Grover Block %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         function Generate_Grover_Block(obj,Question)
             obj.Generate_DxD_Oracle(Question);
             obj.Generate_J();
             obj.Generate_All_Hadamard_Gates();
         end
         
         function Grover_Diffusion_Block(obj)
             itterations =pi/4*sqrt(2^obj.N);
             for idx = 1:itterations
                obj.Apply_Oracle();
                obj.Apply_All_Hadamard();
                obj.Apply_J();
                obj.Apply_All_Hadamard();
             end
         end
         
         function Grovers_Algorithm(obj)
             
             obj.Apply_All_Hadamard();
             obj.Grover_Diffusion_Block();
         end

              %%%%%%%%%%%%%%%%% Display Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       
       function Bar_Chart_Results(obj)
          
          bar(obj.Measurement_Array)
          title('Measurement Results')
           xlabel('States')
           ylabel('Measurements')
       end

         
    end
end

