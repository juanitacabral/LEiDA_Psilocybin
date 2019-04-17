function plot_nodes_in_cortex (V)

% PLOT NODES
hold on

% % PLOT CORTEX
% 
cortex.path='MNI152_T1_2mm_brain_mask.nii';
cortex.pial=mapPial(cortex.path);
cortex.color=[0.9 0.9 0.9];
cortex.transparency=0.1; % To view only opaque cortex =1;
cortex.val=0.2;
redux=1;
sregion=smooth3(cortex.pial);
psregion=patch(isosurface(sregion,cortex.val,'verbose'), 'FaceColor', cortex.color, 'EdgeColor', 'none');
reducepatch(psregion,redux,'verbose');
isonormals(sregion,psregion);
set(psregion,'FaceAlpha', cortex.transparency); %transparency

% Max of V or -V is 1

V=V/max(abs(V));
V=round(V*10)/10;  % round small values close to zero

% center origin
ori=[65 45.5 35];

load aal_cog.txt aal_cog
scale=5.5;
MNI_coord=scale*(aal_cog/10);
clear aal_cog


a=2.5;
[x,y,z] = sphere;
x=a*x;
y=a*y;
z=a*z;

for n=1:length(V)
    if V(n)>0
        surf(x+MNI_coord(n,2)+ori(1), y+MNI_coord(n,1)+ori(2),z+MNI_coord(n,3)+ori(3),'FaceColor',[1 V(n) 0],'EdgeColor','none','FaceAlpha',0.7);
    elseif V(n)<=0
        surf(x+MNI_coord(n,2)+ori(1), y+MNI_coord(n,1)+ori(2),z+MNI_coord(n,3)+ori(3),'FaceColor',[0 -V(n) 1],'EdgeColor','none','FaceAlpha',0.7);
    end
end
  
n_strong=find(V>0.3);
if numel(n_strong)>1
    u=1;
    
    for a=1:numel(n_strong)
        n=n_strong(a);
        for b=1:a
            p=n_strong(b);
            c1=[MNI_coord(n,2)+ori(1) MNI_coord(n,1)+ori(2) MNI_coord(n,3)+ori(3)];
            c2=[MNI_coord(p,2)+ori(1) MNI_coord(p,1)+ori(2) MNI_coord(p,3)+ori(3)];
            plot3([c1(1) c2(1)],[c1(2) c2(2)],[c1(3) c2(3)],'Color','r'); %cmap(IDX(t),:));
            u=u+1;
        end
    end
end
    
axis off;
axis equal


% -------------------------------------------------------
% Setting image properties - light, material, angle
% -------------------------------------------------------
set(gcf,'Renderer', 'OpenGL') % USE UNDER LINUX FOR TRANSPARENCY 
view(3); axis off;
daspect([1 1 1]);
pbaspect([1 1 1]);
set(gca,'CameraViewAngle', 6);
set(gca, 'Projection', 'orthographic')
set(gca, 'CameraTarget', [51 68 90])


%view([-90 60])
%view([90 -90]) % ventral
view([-90 90]) % top
%view([0 0]) % R sideways
% view([-180 0]) % L sideways
% view([45 20]) % perspective
%view([90 0]) % front

material dull; lighting phong;
camlight;
rotate3d;

end

function pial=mapPial(region)

VG=spm_vol(region(1,:));
pial=zeros(VG.dim(1:3)); 
for i=1:VG.dim(3)
  pial(:,:,i) = spm_slice_vol(VG,spm_matrix([0 0 i]),VG.dim(1:2),1);
end

end