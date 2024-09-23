%%Attempt 3 - using pde toolbox
%%Ashwin Krishnan

clf
model = createpde(1); %%Creates a model with 1 pde
gmCylinder = multicylinder(1,4, "Zoffset",2);%%creates cylinder object with arbitrary offset 
gmBox = multicuboid(7,7,7);%%creates cube object
gmTotal = addVoid(gmBox,gmCylinder);%%cylindical hole created in cube


InternalFaces = [7, 8, 9]%cylinder faces
ExternalFaces = 1:1:6%cube faces

hgold = -315e3;%heatflux of electrode
Telectrode = 288;%temperature of electrode

model.Geometry = gmTotal;%add geometry to the pde model

%display the geometry
pdegplot(gmBox,"FaceLabels","on","FaceAlpha",0.25)
hold on
pdegplot(gmCylinder,"CellLabels","on","FaceAlpha",0.25)


pdegplot(model,"FaceLabels","on", "FaceAlpha", 0.5); % "FaceLabels" for 3-D
mesh = generateMesh(model);
applyBoundaryCondition(model, "dirichlet","Face", ExternalFaces, "u", 300);%%Apply dirichet boundary conditions to walls of cube
applyBoundaryCondition(model, "neuman","Face", InternalFaces, "q", hgold*Telectrode,"g", Telectrode);%%Apply neuman boundary conditions to walls of electrode

setInitialConditions(model,310);%%Brain is initially 310K

qm = 49937; %metabolic heat generation in W/m3
pt = 1046; %tissue density
ct = 3630; %tissue specific heat
cb = 3617; %blood specific heat 
pb = 1050; %blood density
wb = 0.02; %perfusion rate (1/s)
Ta = 310; %Arterial blood temperature
kb = 0.51;%brain thermal conductivity


timestep = 5;

%%coefficients for pde solver
f = wb .* pb .* cb .* Ta + qm;
c = kb;
d = pt .* ct;
a = wb .* pb .* cb;
specifyCoefficients(model,"m",0,...
                          "d",d,...
                          "c",c,...
                          "a",a,...
                          "f",f);

tlist = linspace(0,timestep,100);%%time with a timestep interval
results = solvepde(model,tlist)

figure
for i = 0:0.1:10
    [X,Z] = meshgrid(-10:0.1:10,-10:0.1:10);
    Y = i .* ones(size(X));
    V = interpolateSolution(results,X,Y,Z,length(tlist));
    V = reshape(V,size(X));
    subplot(2,1,1)
    contour(X,Z,V);
    axis equal
    title("Contour Plot on Tilted Plane")
    xlabel("z")
    ylabel("y")
    colorbar
    clim([270,350]);
    subplot(2,1,2)
    surf(X,Z,V,"LineStyle","none");
    axis equal
    view(0,90)
    title("Colored Plot on Tilted Plane")
    xlabel("z")
    ylabel("y")
    colorbar
    clim([270,350]);
    drawnow limitrate
    clf
    hold on
end

[X,Z] = meshgrid(-10:0.1:10,-10:0.1:10);
Y = zeros(size(X));
V = interpolateSolution(results,X,Y,Z,length(tlist));
V = reshape(V,size(X));
subplot(2,1,1)
contour(X,Z,V);
axis equal
title("Contour Plot on Tilted Plane")
xlabel("z")
ylabel("y")
colorbar
clim([270,350]);
subplot(2,1,2)
surf(X,Z,V,"LineStyle","none");
axis equal
view(0,90)
title("Colored Plot on Tilted Plane")
xlabel("z")
ylabel("y")
colorbar
clim([270,350]);




