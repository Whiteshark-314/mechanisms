classdef Simple_6R
    properties
        a
        p
        e
        F
        Fx
        Fy
        thetas
        phis
        alpha
        Ex
        Ey
        O
        thetas_combi
        traj_x
        traj_y
        traj_alpha
        traj_thetas_1
        traj_thetas_2
        traj_thetas_3
        r
        ws_x
        ws_y
        
    end
    
    methods
        function obj = Simple_6R(active_links,passive_links,end_effector_sideLengths,fixed_coordinates)
            obj.a=active_links;
            obj.p=passive_links;
            obj.e=end_effector_sideLengths;
            obj.F=fixed_coordinates;
            obj.Fx=obj.F(:,1);
            obj.Fy=obj.F(:,2);
            if 2*(mean(obj.a)+mean(obj.p)+mean(obj.e)/sqrt(3))/sqrt(3) < 1.1*(obj.Fx(2)-obj.Fx(1))
                fprintf("Consider larger values for a,p,e or smaller value for F")
            end
        end
        
        function obj = fk(obj,actuation_angles)
            obj.thetas=actuation_angles;
            
            zeta=pi-acos((obj.e(1)^2+obj.e(2)^2-obj.e(3)^2)/(2*obj.e(1)*obj.e(2)));
            
            H=(obj.Fx(2)-obj.Fx(1)-obj.a(1)*cos(obj.thetas(1))-obj.p(1)*cos(obj.thetas(2))-obj.e(1)*cos(obj.thetas(3)))/obj.p(2);
            I=(obj.Fy(2)-obj.Fy(1)-obj.a(1)*sin(obj.thetas(1))-obj.p(1)*sin(obj.thetas(2))-obj.e(1)*sin(obj.thetas(3)))/obj.p(2);
            J=obj.a(2)/obj.p(2);
            
            c0=H^2+I^2+J^2-2*H*J-1;
            c1=4*I*J;
            c2=H^2+I^2+J^2+2*H*J-1;
            T=roots([c0 c1 c2]);
            T=T(imag(T)==0);
            obj.phis(2,:)=2*atan(T);
            obj.phis(1,:)=atan2(-(I+J*sin(obj.phis(2,:))),-(H+J*cos(obj.phis(2,:))));
            
            obj.Ex(1,:)=obj.Fx(1)+obj.a(1)*cos(obj.thetas(1))+obj.p(1)*cos(obj.thetas(2));
            obj.Ex(2,:)=obj.Fx(2)+obj.a(1)*cos(obj.thetas(1))+obj.p(1)*cos(obj.thetas(2))+obj.e(1)*cos(obj.thetas(3));
            obj.Ex(3,:)=obj.Fx(3)+obj.a(1)*cos(obj.thetas(1))+obj.p(1)*cos(obj.thetas(2))+obj.e(3)*cos(obj.thetas(3)+zeta);
            obj.Ex=obj.Ex';
            obj.Ey(1,:)=obj.Fy(1)+obj.a(1)*sin(obj.thetas(1))+obj.p(1)*sin(obj.thetas(2));
            obj.Ey(2,:)=obj.Fy(2)+obj.a(1)*sin(obj.thetas(1))+obj.p(1)*sin(obj.thetas(2))+obj.e(1)*sin(obj.thetas(3));
            obj.Ey(3,:)=obj.Fy(3)+obj.a(1)*sin(obj.thetas(1))+obj.p(1)*sin(obj.thetas(2))+obj.e(3)*sin(obj.thetas(3)+zeta);
            obj.Ey=obj.Ey';
            obj.O(:,1)=mean(obj.Ex,2);
            obj.O(:,2)=mean(obj.Ey,2);
        end
    end
end