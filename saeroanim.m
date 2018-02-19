function [sys,x0,str,ts] = saeroanim(t,x,u,flag,Config)

switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
     [sys,x0,str,ts]=mdlInitializeSizes(Config);

  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2,
     sys = [];
     if Config.Animenable
       mdlUpdate(t,x,u,Config);
     end;
end

function [sys,x0,str,ts]=mdlInitializeSizes(Config)

%
% Set Sizes Matrix
%
sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 7;
sizes.NumInputs      = 9;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialise the initial conditions
%
x0  = [];

%
% str is always an empty matrix
%
str = [];

%
% initialise the array of sample times
%
%ts  = [-2 0]; % variable sample time
ts  = [Config.update 0]; 
if ~Config.Animenable
   return
end

%
% Initialise Figure Window
% 
     h_anim=figure;
     set(h_anim,'Tag','6DOF anim');


%
% Initialize Axes
%
    handle.axes(1)=axes;
    grid on;
hold on;
xlabel('x')
ylabel('y')
zlabel('z')

[v,f,Nant]=obj3d;
 A0=[1 0 0
     0 1 0
     0 0  1]
v_=(A0*v')';

%
% Initialize Trajectory 
%
   handle.ln1= line(v_(Nant(1),1),v_(Nant(1),2),v_(Nant(1),3),'Color','b');
   handle.ln2= line(v_(Nant(2),1),v_(Nant(2),2),v_(Nant(2),3),'Color','r');
    
%
% Draw in sun Position
%   
Rs=A0*evalin('base','Rsun')*4;
[xs,ys,zs]=sphere(8);
Rds=0;
xs=xs*Rds;
ys=ys*Rds;
zs=zs*Rds;
handle.target=surf(xs+Rs(1)*ones(size(xs)),ys+Rs(2)*ones(size(xs)),zs+Rs(3)*ones(size(xs)),'Cdata',0.65*ones(size(xs)),'FaceLighting','none');
set(handle.target,'facecolor',[1 1 0],'edgecolor',[0 0 0]);
   shading interp; 
   lighting gouraud;
   camlight infinite;     
   handle.sun=Rs;
   
%
% Draw in Missile Shape
%
handle.craft =patch('Vertices',v_,'Faces',f,'FaceVertexCData',([1 1 1]'*0.7*ones(1,length(v)))','FaceColor','flat');
% 
 	handle.vert =v;
    handle.Nant=Nant;
    
    set(handle.craft,'erasemode','nor','edgecolor',[0 0 0],'clipping','off');
    
%
% Set Handles of graphics in Figure UserData
%   
 handle.A0=A0;
set(h_anim,'userdata',handle);
view(3); 
axis square
xlim([-4 4])
ylim([-4 4])
zlim([-4 4]);

lt=light('Position',Rs,'Style','infinite','Color','w');

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function mdlUpdate(t,x,u,Config,count)

%
% Obtain Handles of Figure Objects
%
    handle = get(findobj('type','figure','Tag','6DOF anim'),'userdata');

    if isempty(findobj('type','figure','Tag','6DOF anim'))
     %figure has been manually closed
     return
    end

%
% Form Transformation Matrix
%
      
%
% Update Craft Object 
%
   v = handle.vert;
   Nant=handle.Nant;
   A=handle.A0*reshape(u,3,3)';
v_=(A*v')';
set(handle.craft,'vertices',v_);
lx= get(handle.ln1,'Xdata');  
ly= get(handle.ln1,'Ydata');
lz= get(handle.ln1,'Zdata');
lx=[lx v_(Nant(1),1)];
ly=[ly v_(Nant(1),2)];
lz=[lz v_(Nant(1),3)];
set(handle.ln1,'Xdata',lx,'Ydata',ly,'Zdata',lz);	


lx= get(handle.ln2,'Xdata');  
ly= get(handle.ln2,'Ydata');
lz= get(handle.ln2,'Zdata');
lx=[lx v_(Nant(2),1)];
ly=[ly v_(Nant(2),2)];
lz=[lz v_(Nant(2),3)];
set(handle.ln2,'Xdata',lx,'Ydata',ly,'Zdata',lz);	

%
% Force MATLAB to Update Drawing
%
   drawnow
sys=v_(Nant(1),1:3);
% end mdlUpdate

