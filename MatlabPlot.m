Data = fopen('Output.txt')
s = textscan(Data, '%d %d %f','headerlines',1)
fclose(Data);

X=reshape(s{1},21,201);
size(X)
Y=meshgrid(0:200,0:20);
size(Y)
Z=reshape(s{3},21,201);
surf(X,Y,Z)
