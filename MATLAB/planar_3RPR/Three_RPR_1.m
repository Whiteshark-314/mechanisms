classdef Three_RPR_1
    properties
        rho 
        e
        F
        Fx
        Fy
        thetas
        alpha
        Ex
        Ey
        O
        rho_combi
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
        function obj = Three_RPR_1(end_effector_sideLengths,fixed_coordinates)
            obj.e=end_effector_sideLengths;
            obj.F=fixed_coordinates;
            obj.Fx=obj.F(:,1);
            obj.Fy=obj.F(:,2);
        end
        
        function obj = fk(obj,actuation_distances)
            obj=ik(obj,[sum(obj.Fx)/3,sum(obj.Fy)/3],0);
            min_rho=min(min(abs(obj.rho)));
            obj.Ex=[];
            obj.Ey=[];
            obj.O=[];
            obj.rho=actuation_distances;
            
            zeta=pi-acos((obj.e(1)^2+obj.e(2)^2-obj.e(3)^2)/(2*obj.e(1)*obj.e(2)));
            
            Xx_1=(obj.Fx(2)-obj.Fx(1))/obj.rho(1);
            Yy_1=(obj.Fy(2)-obj.Fy(1))/obj.rho(1);
            pp_1=obj.rho(2)/obj.rho(1);
            ep_1=-obj.e(1)/obj.rho(1);

            Xx_3=(obj.Fx(2)-obj.Fx(3))/obj.rho(3);
            Yy_3=(obj.Fy(2)-obj.Fy(3))/obj.rho(3);
            pp_3=obj.rho(2)/obj.rho(3);
            ep_3=obj.e(2)/obj.rho(3);

            a11=2*Xx_1*pp_1+2*pp_1*ep_1;
            a12=0;
            a13=2*Xx_1*pp_1-2*pp_1*ep_1;
            b11=2*Yy_1*pp_1;
            b12=4*pp_1*ep_1;
            b13=2*Yy_1*pp_1;
            c11=Xx_1^2+Yy_1^2+pp_1^2+ep_1^2-1+2*Xx_1*ep_1;
            c12=4*Yy_1*ep_1;
            c13=Xx_1^2+Yy_1^2+pp_1^2+ep_1^2-1-2*Xx_1*ep_1;

            a31=2*Xx_3*pp_3+2*pp_3*ep_3*cos(zeta);
            a32=-4*pp_3*ep_3*sin(zeta);
            a33=2*Xx_3*pp_3-2*pp_3*ep_3*cos(zeta);
            b31=2*Yy_3*pp_3+2*pp_3*ep_3*sin(zeta);
            b32=4*pp_3*ep_3*cos(zeta);
            b33=2*Yy_3*pp_3-2*pp_3*ep_3*sin(zeta);
            c31=Xx_3^2+Yy_3^2+pp_3^2+ep_3^2-1+2*Xx_3*ep_3*cos(zeta)+2*Yy_3*ep_3*sin(zeta);
            c32=-4*Xx_3*ep_3*sin(zeta)+4*Yy_3*ep_3*cos(zeta);
            c33=Xx_3^2+Yy_3^2+pp_3^2+ep_3^2-1-2*Xx_3*ep_3*cos(zeta)-2*Yy_3*ep_3*sin(zeta);

            p0=((a13*c33 - a33*c13)^2 - (a13*b33 - a33*b13)^2 + (b13*c33 - b33*c13)^2);
            p1=(2*(a13*c33 - a33*c13)*(a12*c33 + a13*c32 - a32*c13 - a33*c12)...
                - 2*(a13*b33 - a33*b13)*(a12*b33 + a13*b32 - a32*b13 - a33*b12)...
                + 2*(b13*c33 - b33*c13)*(b12*c33 + b13*c32 - b32*c13 - b33*c12));
            p2=(2*(a13*c33 - a33*c13)*(a11*c33 + a12*c32 + a13*c31 - a31*c13- a32*c12 - a33*c11)...
                - 2*(a13*b33 - a33*b13)*(a11*b33 + a12*b32 + a13*b31 - a31*b13 - a32*b12 - a33*b11)...
                + 2*(b13*c33 - b33*c13)*(b11*c33 + b12*c32 + b13*c31 - b31*c13 - b32*c12 - b33*c11)...
                - (a12*b33 + a13*b32 - a32*b13 - a33*b12)^2 + (a12*c33 + a13*c32 - a32*c13 - a33*c12)^2 ...
                + (b12*c33 + b13*c32 - b32*c13 - b33*c12)^2);
            p3=(2*(a12*c33 + a13*c32 - a32*c13 - a33*c12)*(a11*c33 + a12*c32+ a13*c31 - a31*c13 - a32*c12 - a33*c11)...
                - 2*(a12*b33 + a13*b32 - a32*b13 - a33*b12)*(a11*b33 + a12*b32 + a13*b31 - a31*b13 - a32*b12 - a33*b11)...
                + 2*(b12*c33 + b13*c32 - b32*c13 - b33*c12)*(b11*c33 + b12*c32 + b13*c31 - b31*c13 - b32*c12 - b33*c11)...
                - 2*(a13*b33 - a33*b13)*(a11*b32 + a12*b31 - a31*b12 - a32*b11)...
                + 2*(a13*c33 - a33*c13)*(a11*c32 + a12*c31 - a31*c12 - a32*c11)...
                + 2*(b13*c33 - b33*c13)*(b11*c32 + b12*c31 - b31*c12 - b32*c11));
            p4=(2*(a11*c32 + a12*c31 - a31*c12 - a32*c11)*(a12*c33 + a13*c32 - a32*c13 - a33*c12)...
                - 2*(a11*b32 + a12*b31 - a31*b12 - a32*b11)*(a12*b33 + a13*b32 - a32*b13 - a33*b12)...
                + 2*(b11*c32 + b12*c31 - b31*c12 - b32*c11)*(b12*c33 + b13*c32 - b32*c13 - b33*c12)...
                - (a11*b33 + a12*b32 + a13*b31 - a31*b13 - a32*b12 - a33*b11)^2 ...
                + (a11*c33 + a12*c32 + a13*c31 - a31*c13 - a32*c12 - a33*c11)^2 ...
                + (b11*c33 + b12*c32 + b13*c31 - b31*c13 - b32*c12 - b33*c11)^2 ...
                - 2*(a11*b31 - a31*b11)*(a13*b33 - a33*b13) + 2*(a11*c31 - a31*c11)*(a13*c33 - a33*c13)...
                + 2*(b11*c31 - b31*c11)*(b13*c33 - b33*c13));
            p5=(2*(a11*c32 + a12*c31 - a31*c12 - a32*c11)*(a11*c33 + a12*c32 + a13*c31 - a31*c13 - a32*c12 - a33*c11)...
                - 2*(a11*b32 + a12*b31 - a31*b12 - a32*b11)*(a11*b33 + a12*b32 + a13*b31 - a31*b13 - a32*b12 - a33*b11)...
                + 2*(b11*c32 + b12*c31 - b31*c12 - b32*c11)*(b11*c33 + b12*c32 + b13*c31 - b31*c13 - b32*c12 - b33*c11)...
                - 2*(a11*b31 - a31*b11)*(a12*b33 + a13*b32 - a32*b13 - a33*b12)...
                + 2*(a11*c31 - a31*c11)*(a12*c33 + a13*c32 - a32*c13 - a33*c12)...
                + 2*(b11*c31 - b31*c11)*(b12*c33 + b13*c32 - b32*c13 - b33*c12));
            p6=(2*(a11*c31 - a31*c11)*(a11*c33 + a12*c32 + a13*c31 - a31*c13 - a32*c12 - a33*c11)...
                - 2*(a11*b31 - a31*b11)*(a11*b33 + a12*b32 + a13*b31 - a31*b13 - a32*b12 - a33*b11)...
                + 2*(b11*c31 - b31*c11)*(b11*c33 + b12*c32 + b13*c31 - b31*c13 - b32*c12 - b33*c11)...
                - (a11*b32 + a12*b31 - a31*b12 - a32*b11)^2 + (a11*c32 + a12*c31 - a31*c12 - a32*c11)^2 ...
                + (b11*c32 + b12*c31 - b31*c12 - b32*c11)^2);
            p7=(2*(a11*c31 - a31*c11)*(a11*c32 + a12*c31 - a31*c12 - a32*c11)...
                - 2*(a11*b31 - a31*b11)*(a11*b32 + a12*b31 - a31*b12 - a32*b11)...
                + 2*(b11*c31 - b31*c11)*(b11*c32 + b12*c31 - b31*c12 - b32*c11));
            p8=(a11*c31 - a31*c11)^2 - (a11*b31 - a31*b11)^2 + (b11*c31 - b31*c11)^2;
            T=roots([p0,p1,p2,p3,p4,p5,p6,p7,p8]);
            T=round(T,6);
            T=T(imag(T)==0);
            if ~(mean(actuation_distances)>(min_rho)) && isempty(T)
                fprintf("Sum of actuation distances is less than recommended\n")
                fprintf("Consider sum of actuation distances of atleast %d\n",min_rho*3)
                return
            end
            obj.alpha=2*atan(T)';
            obj.thetas=zeros(3,length(obj.alpha));

            A1=2*Xx_1*pp_1+2*pp_1*ep_1*cos(obj.alpha);
            B1=2*Yy_1*pp_1+2*pp_1*ep_1*sin(obj.alpha);
            C1=Xx_1^2+Yy_1^2+pp_1^2+ep_1^2+2*ep_1*(Xx_1*cos(obj.alpha)+Yy_1*sin(obj.alpha))-1;
            A3=2*Xx_3*pp_3+2*pp_3*ep_3*cos(obj.alpha+zeta);
            B3=2*Yy_3*pp_3+2*pp_3*ep_3*sin(obj.alpha+zeta);
            C3=Xx_3^2+Yy_3^2+pp_3^2+ep_3^2+2*ep_3*(Xx_3*cos(obj.alpha+zeta)+Yy_3*sin(obj.alpha+zeta))-1;

            CP_2=(B3.*C1-B1.*C3)./(A3.*B1-A1.*B3);
            SP_2=-(A3.*C1-A1.*C3)./(A3.*B1-A1.*B3);
            obj.thetas(2,:)=atan2(SP_2,CP_2);

            CP_1=Xx_1+pp_1*cos(obj.thetas(2,:))+ep_1*cos(obj.alpha);
            SP_1=Yy_1+pp_1*sin(obj.thetas(2,:))+ep_1*sin(obj.alpha);
            obj.thetas(1,:)=atan2(SP_1,CP_1);

            CP_3=Xx_3+pp_3*cos(obj.thetas(2,:))+ep_3*cos(obj.alpha+zeta);
            SP_3=Yy_3+pp_3*sin(obj.thetas(2,:))+ep_3*sin(obj.alpha+zeta);
            obj.thetas(3,:)=atan2(SP_3,CP_3);

            obj.Ex(1,:)=obj.Fx(1)+obj.rho(1)*cos(obj.thetas(1,:));
            obj.Ex(2,:)=obj.Fx(2)+obj.rho(2)*cos(obj.thetas(2,:));
            obj.Ex(3,:)=obj.Fx(3)+obj.rho(3)*cos(obj.thetas(3,:));
            obj.Ex=obj.Ex';
            obj.Ey(1,:)=obj.Fy(1)+obj.rho(1)*sin(obj.thetas(1,:));
            obj.Ey(2,:)=obj.Fy(2)+obj.rho(2)*sin(obj.thetas(2,:));
            obj.Ey(3,:)=obj.Fy(3)+obj.rho(3)*sin(obj.thetas(3,:));
            obj.Ey=obj.Ey';
            obj.O(:,1)=mean(obj.Ex,2);
            obj.O(:,2)=mean(obj.Ey,2);
        end
        
        function obj = ik(obj,O,alpha)
            iota=acos((obj.e(1)^2-obj.e(2)^2+obj.e(3)^2)/(2*obj.e(1)*obj.e(3)));
            A_1=[1,1,1;-1,1,0;-1,0,1];
            A_2=[1,1,1;-1,1,0;-1,0,1];
            b_1=[3*O(1);obj.e(1)*cos(alpha);obj.e(3)*cos(iota+alpha)];
            b_2=[3*O(2);obj.e(1)*sin(alpha);obj.e(3)*sin(iota+alpha)];
            ex=A_1\b_1;
            ey=A_2\b_2;
            obj.Ex=ex';
            obj.Ey=ey';
            eFx=(ex-obj.Fx)';
            eFy=(ey-obj.Fy)';

            A=eFx.^2+eFy.^2;

            T1=[sqrt(A(1));-sqrt(A(1))];
            T2=[sqrt(A(2));-sqrt(A(2))];
            T3=[sqrt(A(3));-sqrt(A(3))];
            T1=T1(imag(T1)==0);
            T2=T2(imag(T2)==0);
            T3=T3(imag(T3)==0);
            if isempty(T1)||isempty(T2)||isempty(T3)
                obj.rho=[];
            else
                obj.rho=[T1';T2';T3'];
                obj.rho=sort(obj.rho,2);
                unique_list=unique(nchoosek([1,2,1,2,1,2],3),'rows');
                obj.rho_combi=[obj.rho(1,unique_list(:,1))', obj.rho(2,unique_list(:,2))', obj.rho(3,unique_list(:,3))'];
            end
        end
    end
end