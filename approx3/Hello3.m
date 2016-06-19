%% Hello World
% Olivier Sete, June 2016

%%
% (Chebfun3 Example: Fun/HelloWorld3.m)

%%
% Here is the matrix from the Chebfun2 example 'Hello World'
% (see http://www.chebfun.org/examples/fun/HelloWorld.html).
A=zeros(15,40);
A(2:9,2:3)=1; A(5:6,4:5)=1;A(2:9,6:7)=1; A(3:10,10:11)=1;
A(3:4,10:15)=1; A(6:7,10:15)=1; A(9:10,10:15)=1; A(4:11,18:19)=1;
A(10:11,18:24)=1; A(5:12,26:27)=1; A(11:12,26:31)=1;
A(6:13,34:35)=1; A(6:13,38:39)=1; A(6:7,36:37)=1; A(12:13,36:37)=1;

%% 
% We pad it with zeros to size 40x40 and add a third dimension:

A = [zeros(14,40); A; zeros(11,40)];
A = fliplr(flipud(A));

B = zeros(40,40,40);
for k = 18:21
    B(k,:,:) = A;
end

%%
% B is now a 40x40x40 3D tensor whose entries are all zeros and ones.

%% Constructing a chebfun3 from discrete data
% A chebfun3 is normally constructed from a function handle for a 
% function $f(x,y,z)$.
% It is also possible to do it from discrete data, and e can specify that
% this it to be interpreted as living on an equispaced grid:

f = chebfun3(B, 'equi');

%% 
% Let's plot the result:

f = permute(f, [1, 3, 2]);
isosurface(f, 0.5)
view([-2.5, -1, 0.4]), camlight, axis off, axis equal
