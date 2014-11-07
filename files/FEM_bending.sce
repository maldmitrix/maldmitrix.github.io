m = 100;                    // число конечных элементов
P = -100;                   // нагрузка

// Массив из m+1 равномерно распределенных точек на интервале [0,1]
nodeCoordinates = linspace(0, 1, m+1)';
L = max(nodeCoordinates);   // Наибольшая координата узлов, т.е. L=1
numberNodes = size(nodeCoordinates,1);
xx = nodeCoordinates;
E = 1; I = 1; EI = E*I;         // Изгибная жесткость балки
GDof = 2*numberNodes;           // Глобальное число степеней свобод
U = zeros(GDof,1);
force = zeros(GDof,1);          // Глобальный столбец сил
stiffness = zeros(GDof,GDof);   // Глобальная матрица жесткости 
displacements = zeros(GDof,1);  // Глобальный столбец прогибов

// Нумерация узлов каждого элемента (для первого индексы: [1, 2])
for i = 1:m
  elementNodes(i,1) = i;
  elementNodes(i,2) = i+1;
end

for e=1:m                       // Пробегаем все элементы
  indice=elementNodes(e,:);     // Номера узлов элемента (2 штуки)
  // Номера степеней свобод, принадлежащих элементу (4 штуки)
  elementDof = [ 2*(indice(1)-1)+1 2*(indice(2)-1) ...
                 2*(indice(2)-1)+1 2*(indice(2)-1)+2];
  h = xx(indice(2)) - xx(indice(1));  // Длина элемента
  // Набиваем локальный столбец сил каждого элемента
  f1 = (P*h/2)*[1  6*h  1  -6*h]';  
  force(elementDof)=force(elementDof)+f1;  
  // Набиваем локальную матрицу жесткости каждого элемента
  k1=(EI/h^3)*[  12    6*h   -12    6*h ;...
                6*h  4*h^2  -6*h  2*h^2 ;...
                -12   -6*h    12   -6*h ;...
                6*h  2*h^2  -6*h  4*h^2 ]';
  stiffness(elementDof,elementDof)=stiffness(elementDof,elementDof)+k1;
end


// Граничные условия: оба конца защемлены
fixedNodeU = [1 2*m+1];
fixedNodeV = [2 2*m+2];
prescribedDof = [fixedNodeU;fixedNodeV];
activeDof = setdiff([1:GDof]',[prescribedDof]);
U = stiffness(activeDof,activeDof) \ force(activeDof);
displacements = zeros(GDof,1);
displacements(activeDof) = U;
U = displacements(1:2:2*numberNodes);
plot(nodeCoordinates,U,'b:');       // синий штрих чере две точки

// Граничные условия: оба конца свободно оперты
fixedNodeU = [1 2*m+1];
fixedNodeV = [];
prescribedDof = [fixedNodeU;fixedNodeV];
activeDof = setdiff([1:GDof]',[prescribedDof]);
U = stiffness(activeDof,activeDof) \ force(activeDof);
displacements = zeros(GDof,1);
displacements(activeDof) = U;
U = displacements(1:2:2*numberNodes);
plot(nodeCoordinates,U,'r--');      // красный штрих

// Граничные условия: левый конец защемлен
// правый конец свободно оперт
fixedNodeU = [1 2*m+1];
fixedNodeV = [2 2*m+1];
prescribedDof = [fixedNodeU;fixedNodeV];
activeDof = setdiff([1:GDof]',[prescribedDof]);
U = stiffness(activeDof,activeDof) \ force(activeDof);
displacements = zeros(GDof,1);
displacements(activeDof) = U;
U = displacements(1:2:2*numberNodes);
plot(nodeCoordinates,U,'k-'); // черная сплошная линия

// Граничные условия: левая половина балки защемлена, 
// правая половина вся свободна
fixedNodeU = [1:m];
fixedNodeV = [1:m];
prescribedDof = [fixedNodeU;fixedNodeV];
activeDof = setdiff([1:GDof]',[prescribedDof]);
U = stiffness(activeDof,activeDof) \ force(activeDof);
displacements = zeros(GDof,1);
displacements(activeDof) = U;
U = displacements(1:2:2*numberNodes);
plot(nodeCoordinates,U,'g-');       // зеленая сплошная линия 