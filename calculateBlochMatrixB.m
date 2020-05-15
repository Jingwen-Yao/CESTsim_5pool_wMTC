function [b, B, Binvb] = ...
    calculateBlochMatrixB(M0List, KxwList, KwxList, T1List, T2List, dList)

R1List = 1./T1List;
R1List(T1List<1e-15) = 0;
R2List = 1./T2List;
R2List(T2List<1e-15) = 0;

b = [0;0;M0List(1)*R1List(1); ...
    0;0;M0List(2)*R1List(2); ...
    0;0;M0List(3)*R1List(3); ...
    0;0;M0List(4)*R1List(4); ...
    0;0;M0List(5)*R1List(5); ...
    M0List(6)*R1List(6)];

Asub = zeros(3,3,5,5);
for i = 1:5
    Asub(:,:,i,i) = [-R2List(i)-KxwList(i)    0   0;...
                     0  -R2List(i)-KxwList(i)     0;...
                     0          0     -R1List(i)-KxwList(i)];
end
for i = 2:5
    Asub(:,:,1,i) = [KxwList(i) 0           0;...
                     0          KxwList(i)  0;...
                     0          0           KxwList(i)];
    Asub(:,:,i,1) = [KwxList(i) 0           0;...
                     0          KwxList(i)  0;...
                     0          0           KwxList(i)];
end
A = [];
for i = 1:5
    Arow = [Asub(:,:,i,1) Asub(:,:,i,2) Asub(:,:,i,3) Asub(:,:,i,4) Asub(:,:,i,5)];
    A = [A;Arow];
end

A = [A [0;0;KxwList(6);0;0;0;0;0;0;0;0;0;0;0;0]];
B = [A;[0 0 KwxList(6) 0 0 0 0 0 0 0 0 0 0 0 0 R1List(6)]];
B(abs(B)<1e-15) = 0;
B(isnan(B)) = 0;

Binvb = B\b;

end