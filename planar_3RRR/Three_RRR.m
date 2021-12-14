classdef Three_RRR
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
        delta
        lambdas
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
        function obj = Three_RRR(active_links,passive_links,...
                                end_effector_sideLengths,fixed_coordinates)
            obj.a=active_links;
            obj.p=passive_links;
            obj.e=end_effector_sideLengths;
            obj.F=fixed_coordinates;
            obj.Fx=obj.F(:,1);
            obj.Fy=obj.F(:,2);
        end

        function obj = fk(obj,actuation_angles,delta,lambdas)
            %obj=ik(obj,[sum(obj.Fx)/3,sum(obj.Fy)/3],0);
            %min_rho=min(min(abs(obj.rho)));
            obj.Ex=[];
            obj.Ey=[];
            obj.O=[];
            if nargin==4
                obj.delta=delta;
                lambdas=reshape(lambdas,1,9,[]);
                N=size(lambdas,3);
            elseif size(obj.a,1)>1
                obj.delta=0;
                lambdas=zeros(1,9);
                N=1;
                P=size(obj.a,1);
            else
                obj.delta=0;
                lambdas=zeros(1,9);
                N=1;
                P=size(actuation_angles,1);
            end
            obj.thetas=actuation_angles;

            zeta=pi-acos((obj.e(:,1).^2+obj.e(:,2).^2-...
                obj.e(:,3).^2)./(2*obj.e(:,1).*obj.e(:,2)));

            del_C1=obj.delta*(cos(lambdas(1,2,:))+cos(lambdas(1,5,:))+...
                cos(lambdas(1,8,:))-cos(lambdas(1,1,:))-...
                cos(lambdas(1,4,:))-cos(lambdas(1,7,:)));
            del_S1=obj.delta*(sin(lambdas(1,2,:))+sin(lambdas(1,5,:))+...
                sin(lambdas(1,8,:))-sin(lambdas(1,1,:))-...
                sin(lambdas(1,4,:))-sin(lambdas(1,7,:)));
            Xx_1=(obj.Fx(2)-obj.Fx(1)+obj.a(2)*cos(obj.thetas(:,2))-...
                obj.a(:,1)*cos(obj.thetas(:,1))+del_C1)./obj.p(:,1);
            Yy_1=(obj.Fy(2)-obj.Fy(1)+obj.a(2)*sin(obj.thetas(:,2))-...
                obj.a(:,1)*sin(obj.thetas(:,1))+del_S1)./obj.p(:,1);
            pp_1=obj.p(:,2)./obj.p(:,1);
            ep_1=-obj.e(:,1)./obj.p(:,1);

            del_C3=obj.delta*(cos(lambdas(1,2,:))+cos(lambdas(1,5,:))+...
                cos(lambdas(1,8,:))-cos(lambdas(1,3,:))-...
                cos(lambdas(1,6,:))-cos(lambdas(1,9,:)));
            del_S3=obj.delta*(sin(lambdas(1,2,:))+sin(lambdas(1,5,:))+...
                sin(lambdas(1,8,:))-sin(lambdas(1,3,:))-...
                sin(lambdas(1,6,:))-sin(lambdas(1,9,:)));
            Xx_3=(obj.Fx(2)-obj.Fx(3)+obj.a(2)*cos(obj.thetas(:,2))-...
                obj.a(:,3)*cos(obj.thetas(:,3))+del_C3)./obj.p(:,3);
            Yy_3=(obj.Fy(2)-obj.Fy(3)+obj.a(2)*sin(obj.thetas(:,2))-...
                obj.a(:,3)*sin(obj.thetas(:,3))+del_S3)./obj.p(:,3);
            pp_3=obj.p(:,2)./obj.p(:,3);
            ep_3=obj.e(:,2)./obj.p(:,3);

            a11=2*Xx_1.*pp_1+2*pp_1.*ep_1;
            a12=0;
            a13=2*Xx_1.*pp_1-2*pp_1.*ep_1;
            b11=2*Yy_1.*pp_1;
            b12=4*pp_1.*ep_1;
            b13=2*Yy_1.*pp_1;
            c11=Xx_1.^2+Yy_1.^2+pp_1.^2+ep_1.^2-1+2*Xx_1.*ep_1;
            c12=4*Yy_1.*ep_1;
            c13=Xx_1.^2+Yy_1.^2+pp_1.^2+ep_1.^2-1-2*Xx_1.*ep_1;

            a31=2*Xx_3.*pp_3+2*pp_3.*ep_3.*cos(zeta);
            a32=-4*pp_3.*ep_3.*sin(zeta);
            a33=2*Xx_3.*pp_3-2*pp_3.*ep_3.*cos(zeta);
            b31=2*Yy_3.*pp_3+2*pp_3.*ep_3.*sin(zeta);
            b32=4*pp_3.*ep_3.*cos(zeta);
            b33=2*Yy_3.*pp_3-2*pp_3.*ep_3.*sin(zeta);
            c31=Xx_3.^2+Yy_3.^2+pp_3.^2+ep_3.^2-1+2*Xx_3.*ep_3.*cos(zeta)+...
                2*Yy_3.*ep_3.*sin(zeta);
            c32=-4*Xx_3.*ep_3.*sin(zeta)+4*Yy_3.*ep_3.*cos(zeta);
            c33=Xx_3.^2+Yy_3.^2+pp_3.^2+ep_3.^2-1-2*Xx_3.*ep_3.*cos(zeta)-...
                2*Yy_3.*ep_3.*sin(zeta);

            p0=((a13.*c33 - a33.*c13).^2-(a13.*b33 - a33.*b13).^2+(b13.*c33 - b33.*c13).^2);
            p1=(2*(a13.*c33 - a33.*c13).*(a12.*c33 + a13.*c32 - a32.*c13 - a33.*c12)...
                - 2*(a13.*b33 - a33.*b13).*(a12.*b33 + a13.*b32 - a32.*b13 - a33.*b12)...
                + 2*(b13.*c33 - b33.*c13).*(b12.*c33 + b13.*c32 - b32.*c13 - b33.*c12));
            p2=(2*(a13.*c33 - a33.*c13).*(a11.*c33 + a12.*c32 + a13.*c31 - a31.*c13- a32.*c12 - a33.*c11)...
                - 2*(a13.*b33 - a33.*b13).*(a11.*b33 + a12.*b32 + a13.*b31 - a31.*b13 - a32.*b12 - a33.*b11)...
                + 2*(b13.*c33 - b33.*c13).*(b11.*c33 + b12.*c32 + b13.*c31 - b31.*c13 - b32.*c12 - b33.*c11)...
                - (a12.*b33 + a13.*b32 - a32.*b13 - a33.*b12).^2 + (a12.*c33 + a13.*c32 - a32.*c13 - a33.*c12).^2 ...
                + (b12.*c33 + b13.*c32 - b32.*c13 - b33.*c12).^2);
            p3=(2*(a12.*c33 + a13.*c32 - a32.*c13 - a33.*c12).*(a11.*c33 + a12.*c32+ a13.*c31 - a31.*c13 - a32.*c12 - a33.*c11)...
                - 2*(a12.*b33 + a13.*b32 - a32.*b13 - a33.*b12).*(a11.*b33 + a12.*b32 + a13.*b31 - a31.*b13 - a32.*b12 - a33.*b11)...
                + 2*(b12.*c33 + b13.*c32 - b32.*c13 - b33.*c12).*(b11.*c33 + b12.*c32 + b13.*c31 - b31.*c13 - b32.*c12 - b33.*c11)...
                - 2*(a13.*b33 - a33.*b13).*(a11.*b32 + a12.*b31 - a31.*b12 - a32.*b11)...
                + 2*(a13.*c33 - a33.*c13).*(a11.*c32 + a12.*c31 - a31.*c12 - a32.*c11)...
                + 2*(b13.*c33 - b33.*c13).*(b11.*c32 + b12.*c31 - b31.*c12 - b32.*c11));
            p4=(2*(a11.*c32 + a12.*c31 - a31.*c12 - a32.*c11).*(a12.*c33 + a13.*c32 - a32.*c13 - a33.*c12)...
                - 2*(a11.*b32 + a12.*b31 - a31.*b12 - a32.*b11).*(a12.*b33 + a13.*b32 - a32.*b13 - a33.*b12)...
                + 2*(b11.*c32 + b12.*c31 - b31.*c12 - b32.*c11).*(b12.*c33 + b13.*c32 - b32.*c13 - b33.*c12)...
                - (a11.*b33 + a12.*b32 + a13.*b31 - a31.*b13 - a32.*b12 - a33.*b11).^2 ...
                + (a11.*c33 + a12.*c32 + a13.*c31 - a31.*c13 - a32.*c12 - a33.*c11).^2 ...
                + (b11.*c33 + b12.*c32 + b13.*c31 - b31.*c13 - b32.*c12 - b33.*c11).^2 ...
                - 2*(a11.*b31 - a31.*b11).*(a13.*b33 - a33.*b13) + 2*(a11.*c31 - a31.*c11).*(a13.*c33 - a33.*c13)...
                + 2*(b11.*c31 - b31.*c11).*(b13.*c33 - b33.*c13));
            p5=(2*(a11.*c32 + a12.*c31 - a31.*c12 - a32.*c11).*(a11.*c33 + a12.*c32 + a13.*c31 - a31.*c13 - a32.*c12 - a33.*c11)...
                - 2*(a11.*b32 + a12.*b31 - a31.*b12 - a32.*b11).*(a11.*b33 + a12.*b32 + a13.*b31 - a31.*b13 - a32.*b12 - a33.*b11)...
                + 2*(b11.*c32 + b12.*c31 - b31.*c12 - b32.*c11).*(b11.*c33 + b12.*c32 + b13.*c31 - b31.*c13 - b32.*c12 - b33.*c11)...
                - 2*(a11.*b31 - a31.*b11).*(a12.*b33 + a13.*b32 - a32.*b13 - a33.*b12)...
                + 2*(a11.*c31 - a31.*c11).*(a12.*c33 + a13.*c32 - a32.*c13 - a33.*c12)...
                + 2*(b11.*c31 - b31.*c11).*(b12.*c33 + b13.*c32 - b32.*c13 - b33.*c12));
            p6=(2*(a11.*c31 - a31.*c11).*(a11.*c33 + a12.*c32 + a13.*c31 - a31.*c13 - a32.*c12 - a33.*c11)...
                - 2*(a11.*b31 - a31.*b11).*(a11.*b33 + a12.*b32 + a13.*b31 - a31.*b13 - a32.*b12 - a33.*b11)...
                + 2*(b11.*c31 - b31.*c11).*(b11.*c33 + b12.*c32 + b13.*c31 - b31.*c13 - b32.*c12 - b33.*c11)...
                - (a11.*b32 + a12.*b31 - a31.*b12 - a32.*b11).^2 + (a11.*c32 + a12.*c31 - a31.*c12 - a32.*c11).^2 ...
                + (b11.*c32 + b12.*c31 - b31.*c12 - b32.*c11).^2);
            p7=(2*(a11.*c31 - a31.*c11).*(a11.*c32 + a12.*c31 - a31.*c12 - a32.*c11)...
                - 2*(a11.*b31 - a31.*b11).*(a11.*b32 + a12.*b31 - a31.*b12 - a32.*b11)...
                + 2*(b11.*c31 - b31.*c11).*(b11.*c32 + b12.*c31 - b31.*c12 - b32.*c11));
            p8=(a11.*c31 - a31.*c11).^2 - (a11.*b31 - a31.*b11).^2 + (b11.*c31 - b31.*c11).^2;
            T=[];
            if N>1
                for i=N:-1:1
                    t(1,:,i)=roots([p0(1,1,i),p1(1,1,i),p2(1,1,i),p3(1,1,i),p4(1,1,i),p5(1,1,i),p6(1,1,i),p7(1,1,i),p8(1,1,i)])';
                    t(1,:,i)=round(t(1,:,i),10);
                    temp=t(1,:,i);
                    if sum(imag(t(1,:,i))==0)==0
                        T(1,:,i)=[NaN, NaN];
                    else
                        temp2=temp(imag(t(1,:,i))==0);
                        T(1,:,i)=temp2(1:2);
                    end
                end
                obj.alpha=sort(2*atan(T));
                obj.phis=zeros(3,length(obj.alpha));
            end
            if P>1
                for i=P:-1:1
                    tt=complex(NaN(1,8));
                    R=roots([p0(i,1),p1(i,1),p2(i,1),p3(i,1),p4(i,1),p5(i,1),p6(i,1),p7(i,1),p8(i,1)])';
                    l=length(R);
                    tt(1:l)=R;
                    t(i,:,1)=tt; 
                    t(i,:,1)=round(t(i,:,1),10);
                    temp=t(i,:,1);
                    if sum(imag(t(i,:,1))==0)==0
                        T(i,:,1)=NaN(1,8);
                    else
                        temp2=temp(imag(t(i,:,1))==0);
                        L=length(temp2);
                        T(i,:,1)=[temp2,NaN(1,8-L)];
                    end
                end
                obj.alpha=sort(2*atan(T),2);
                obj.phis=zeros(length(obj.alpha),8,3);
            end
            if P==1 && N==1
                R=roots([p0, p1, p2, p3, p4, p5, p6, p7, p8])';
                temp2=R(imag(R)==0);
                L=length(temp2);
                T=[temp2, NaN(1,8-L)];
                obj.alpha=sort(2*atan(T),2);
                obj.phis=zeros(1,8,3);
            end
            if isempty(T)
                return
            end

            A1=2*Xx_1.*pp_1+2*pp_1.*ep_1.*cos(obj.alpha);
            B1=2*Yy_1.*pp_1+2*pp_1.*ep_1.*sin(obj.alpha);
            C1=Xx_1.^2+Yy_1.^2+pp_1.^2+ep_1.^2+2*ep_1.*(Xx_1.*cos(obj.alpha)+Yy_1.*sin(obj.alpha))-1;
            A3=2*Xx_3.*pp_3+2*pp_3.*ep_3.*cos(obj.alpha+zeta);
            B3=2*Yy_3.*pp_3+2*pp_3.*ep_3.*sin(obj.alpha+zeta);
            C3=Xx_3.^2+Yy_3.^2+pp_3.^2+ep_3.^2+2*ep_3.*(Xx_3.*cos(obj.alpha+zeta)+Yy_3.*sin(obj.alpha+zeta))-1;

            CP_2=(B3.*C1-B1.*C3)./(A3.*B1-A1.*B3);
            SP_2=-(A3.*C1-A1.*C3)./(A3.*B1-A1.*B3);
            obj.phis(:,:,2)=atan2(SP_2,CP_2);

            CP_1=Xx_1+pp_1.*cos(obj.phis(:,:,2))+ep_1.*cos(obj.alpha);
            SP_1=Yy_1+pp_1.*sin(obj.phis(:,:,2))+ep_1.*sin(obj.alpha);
            obj.phis(:,:,1)=atan2(SP_1,CP_1);

            CP_3=Xx_3+pp_3.*cos(obj.phis(:,:,2))+ep_3.*cos(obj.alpha+zeta);
            SP_3=Yy_3+pp_3.*sin(obj.phis(:,:,2))+ep_3.*sin(obj.alpha+zeta);
            obj.phis(:,:,3)=atan2(SP_3,CP_3);
            if N>1
                obj.Ex(1,:,:)=obj.Fx(1)+obj.a(1)*cos(obj.thetas(1))+obj.p(1)*cos(obj.phis(1,:))+...
                    obj.delta*(cos(lambdas(1,1,:))+cos(lambdas(1,4,:))+cos(lambdas(1,7,:)));
                obj.Ex(2,:,:)=obj.Fx(2)+obj.a(2)*cos(obj.thetas(2))+obj.p(2)*cos(obj.phis(2,:))+...
                    obj.delta*(cos(lambdas(1,2,:))+cos(lambdas(1,5,:))+cos(lambdas(1,8,:)));
                obj.Ex(3,:,:)=obj.Fx(3)+obj.a(3)*cos(obj.thetas(3))+obj.p(3)*cos(obj.phis(3,:))+...
                    obj.delta*(cos(lambdas(1,3,:))+cos(lambdas(1,6,:))+cos(lambdas(1,9,:)));
                obj.Ey(1,:,:)=obj.Fy(1)+obj.a(1)*sin(obj.thetas(1))+obj.p(1)*sin(obj.phis(1,:))+...
                    obj.delta*(sin(lambdas(1,1,:))+sin(lambdas(1,4,:))+sin(lambdas(1,7,:)));
                obj.Ey(2,:,:)=obj.Fy(2)+obj.a(2)*sin(obj.thetas(2))+obj.p(2)*sin(obj.phis(2,:))+...
                    obj.delta*(sin(lambdas(1,2,:))+sin(lambdas(1,5,:))+sin(lambdas(1,8,:)));
                obj.Ey(3,:,:)=obj.Fy(3)+obj.a(3)*sin(obj.thetas(3))+obj.p(3)*sin(obj.phis(3,:))+...
                    obj.delta*(sin(lambdas(1,3,:))+sin(lambdas(1,6,:))+sin(lambdas(1,9,:)));
                obj.O(1,:,:)=mean(obj.Ex,1);
                obj.O(2,:,:)=mean(obj.Ey,1);
            elseif P>1
                obj.Ex(:,:,1)=obj.Fx(1)+obj.a(1)*cos(obj.thetas(1))+obj.p(1)*cos(obj.phis(:,:,1))+...
                    obj.delta*(cos(lambdas(1,1))+cos(lambdas(1,4))+cos(lambdas(1,7)));
                obj.Ex(:,:,2)=obj.Fx(2)+obj.a(2)*cos(obj.thetas(2))+obj.p(2)*cos(obj.phis(:,:,2))+...
                    obj.delta*(cos(lambdas(1,2))+cos(lambdas(1,5))+cos(lambdas(1,8)));
                obj.Ex(:,:,3)=obj.Fx(3)+obj.a(3)*cos(obj.thetas(3))+obj.p(3)*cos(obj.phis(:,:,3))+...
                    obj.delta*(cos(lambdas(1,3))+cos(lambdas(1,6))+cos(lambdas(1,9)));
                obj.Ey(:,:,1)=obj.Fy(1)+obj.a(1)*sin(obj.thetas(1))+obj.p(1)*sin(obj.phis(:,:,1))+...
                    obj.delta*(sin(lambdas(1,1))+sin(lambdas(1,4))+sin(lambdas(1,7)));
                obj.Ey(:,:,2)=obj.Fy(2)+obj.a(2)*sin(obj.thetas(2))+obj.p(2)*sin(obj.phis(:,:,2))+...
                    obj.delta*(sin(lambdas(1,2))+sin(lambdas(1,5))+sin(lambdas(1,8)));
                obj.Ey(:,:,3)=obj.Fy(3)+obj.a(3)*sin(obj.thetas(3))+obj.p(3)*sin(obj.phis(:,:,3))+...
                    obj.delta*(sin(lambdas(1,3))+sin(lambdas(1,6))+sin(lambdas(1,9)));
                obj.O(:,:,1)=mean(obj.Ex,3);
                obj.O(:,:,2)=mean(obj.Ey,3);
            else
                obj.Ex(1,:)=obj.Fx(1)+obj.a(1)*cos(obj.thetas(1))+obj.p(1)*squeeze(cos(obj.phis(:,:,1)))+...
                    obj.delta*(cos(lambdas(1,1))+cos(lambdas(1,4))+cos(lambdas(1,7)));
                obj.Ex(2,:)=obj.Fx(2)+obj.a(2)*cos(obj.thetas(2))+obj.p(2)*squeeze(cos(obj.phis(:,:,2)))+...
                    obj.delta*(cos(lambdas(1,2))+cos(lambdas(1,5))+cos(lambdas(1,8)));
                obj.Ex(3,:)=obj.Fx(3)+obj.a(3)*cos(obj.thetas(3))+obj.p(3)*squeeze(cos(obj.phis(:,:,3)))+...
                    obj.delta*(cos(lambdas(1,3))+cos(lambdas(1,6))+cos(lambdas(1,9)));
                obj.Ey(1,:)=obj.Fy(1)+obj.a(1)*sin(obj.thetas(1))+obj.p(1)*squeeze(sin(obj.phis(:,:,1)))+...
                    obj.delta*(sin(lambdas(1,1))+sin(lambdas(1,4))+sin(lambdas(1,7)));
                obj.Ey(2,:)=obj.Fy(2)+obj.a(2)*sin(obj.thetas(2))+obj.p(2)*squeeze(sin(obj.phis(:,:,2)))+...
                    obj.delta*(sin(lambdas(1,2))+sin(lambdas(1,5))+sin(lambdas(1,8)));
                obj.Ey(3,:)=obj.Fy(3)+obj.a(3)*sin(obj.thetas(3))+obj.p(3)*squeeze(sin(obj.phis(:,:,3)))+...
                    obj.delta*(sin(lambdas(1,3))+sin(lambdas(1,6))+sin(lambdas(1,9)));
                obj.O(1,:)=mean(obj.Ex,1);
                obj.O(2,:)=mean(obj.Ey,1);
            end
        end

        function obj = ik(obj,O,alpha,delta,lambdas)
            if nargin==5
                obj.delta=delta;
                obj.lambdas=lambdas;
            else
                obj.delta=0;
                obj.lambdas=zeros(1,9);
            end
            iota=acos((obj.e(1)^2-obj.e(2)^2+obj.e(3)^2)/(2*obj.e(1)*obj.e(3)));
            A_1=[1,1,1;-1,1,0;-1,0,1];
            A_2=[1,1,1;-1,1,0;-1,0,1];
            b_1=[3*O(1);obj.e(1)*cos(alpha);obj.e(3)*cos(iota+alpha)];
            b_2=[3*O(2);obj.e(1)*sin(alpha);obj.e(3)*sin(iota+alpha)];
            ex=A_1\b_1;
            ey=A_2\b_2;
            obj.Ex=ex';
            obj.Ey=ey';
            L=reshape(obj.lambdas,[3,3]);
            C=cos(L);
            S=sin(L);
            eFx=(ex-obj.Fx)'+C(:,2)-C(:,1);
            eFy=(ey-obj.Fy)'+S(:,2)-S(:,1);

            A=eFx.^2+eFy.^2+obj.a.^2+2*eFx.*obj.a-obj.p.^2;
            B=-4*eFy.*obj.a;
            C=eFx.^2+eFy.^2+obj.a.^2-2*eFx.*obj.a-obj.p.^2;

            T1=roots([A(1) B(1) C(1)]);
            T2=roots([A(2) B(2) C(2)]);
            T3=roots([A(3) B(3) C(3)]);
            T1=T1(imag(T1)==0);
            T2=T2(imag(T2)==0);
            T3=T3(imag(T3)==0);
            if isempty(T1)||isempty(T2)||isempty(T3)
                obj.thetas=[];
            else
                obj.thetas=[2*atan(T1)';2*atan(T2)';2*atan(T3)'];
                obj.phis(1,:)=atan((obj.Ex(1)-obj.a(1)*cos(obj.thetas(1,:)))./(obj.Ey(1)-obj.a(1)*cos(obj.thetas(1,:))));
                obj.phis(2,:)=atan((obj.Ex(2)-obj.a(2)*cos(obj.thetas(2,:)))./(obj.Ey(2)-obj.a(2)*cos(obj.thetas(2,:))));
                obj.phis(3,:)=atan((obj.Ex(3)-obj.a(3)*cos(obj.thetas(3,:)))./(obj.Ey(3)-obj.a(3)*cos(obj.thetas(3,:))));
                obj.thetas(obj.thetas<0)=2*pi+obj.thetas(obj.thetas<0);
                obj.thetas=sort(obj.thetas,2);
                unique_list=unique(nchoosek([1,2,1,2,1,2],3),'rows');
                obj.thetas_combi=[obj.thetas(1,unique_list(:,1))', obj.thetas(2,unique_list(:,2))', obj.thetas(3,unique_list(:,3))'];
            end
        end

        function obj=trajectory(obj,x_values,y_values,alpha_values)
            obj.traj_x=x_values;
            obj.traj_y=y_values;
            obj.traj_alpha=alpha_values;
            l=length(x_values);
            obj.traj_thetas_1=zeros(l,8);
            obj.traj_thetas_2=zeros(l,8);
            obj.traj_thetas_3=zeros(l,8);
            t_1=zeros(2,2);
            t_2=zeros(2,2);
            t_3=zeros(2,2);
            for i=1:1:l
                x=obj.traj_x(i);
                y=obj.traj_y(i);
                Alpha=obj.traj_alpha(i);
                obj=ik(obj,[x,y],Alpha);
                if isempty(obj.thetas)
                    continue
                end
                t_1(2,:)=obj.thetas(1,:);
                t_2(2,:)=obj.thetas(2,:);
                t_3(2,:)=obj.thetas(3,:);
                obj.traj_thetas_1(i,:)=obj.thetas_combi(:,1)';
                obj.traj_thetas_2(i,:)=obj.thetas_combi(:,2)';
                obj.traj_thetas_3(i,:)=obj.thetas_combi(:,3)';
                if i>=2
                    t_1(1,:)=obj.traj_thetas_1(i-1,4:5);
                    t_2(1,:)=obj.traj_thetas_2(i-1,2:3);
                    t_3(1,:)=obj.traj_thetas_3(i-1,1:2);
                    if (abs(diff([t_1(1,1),t_1(2,2)]))<0.2)  ||  (abs(diff([t_1(1,2),t_1(2,1)]))<0.2)
                        obj.traj_thetas_1(i,:)=circshift(obj.traj_thetas_1(i,:),4);
                    end
                    if (abs(diff([t_2(1,1),t_2(2,2)]))<0.2)  ||  (abs(diff([t_2(1,2),t_2(2,1)]))<0.2)
                        obj.traj_thetas_2(i,:)=circshift(obj.traj_thetas_2(i+1,:),2);
                    end
                    if (abs(diff([t_3(1,1),t_3(2,2)]))<0.2)  ||  (abs(diff([t_3(1,2),t_3(2,1)]))<0.2)
                        obj.traj_thetas_3(i,:)=circshift(obj.traj_thetas_3(i,:),1);
                    end
                else
                    t_1(1,:)=obj.thetas(1,:);
                    t_2(1,:)=obj.thetas(2,:);
                    t_3(1,:)=obj.thetas(3,:);
                end
            end
        end

        function obj = traj_repeat(obj,theta_1_vals,theta_2_vals,theta_3_vals,Alpha)
            obj.traj_thetas_1=theta_1_vals;
            obj.traj_thetas_2=theta_2_vals;
            obj.traj_thetas_3=theta_3_vals;
            l=length(theta_1_vals);
            obj.traj_x=zeros(l,1);
            obj.traj_y=zeros(l,1);
            obj.traj_alpha=zeros(l,1);
            for i=1:l
                obj=fk(obj,[obj.traj_thetas_1(i), obj.traj_thetas_2(i), obj.traj_thetas_3(i)]);
                [~,idx]=min(abs(obj.alpha-Alpha(i)));
                obj.traj_x(i)=obj.O(idx,1);
                obj.traj_y(i)=obj.O(idx,2);
                obj.traj_alpha(i)=obj.alpha(idx);
                obj.Ex=[];
                obj.Ey=[];
                obj.O=[];
            end
        end

        function obj = workspace(obj)
            [k,l]=circcirc(0,0,obj.e(3),obj.e(1),0,obj.e(2));
            e_xy=[0,0;obj.e(1),0;k(1),abs(l(1))];
            O_xy=mean(e_xy);
            obj.r=[obj.a(1)+obj.p(1)+sqrt(sum((e_xy(1,:)-O_xy).^2)) obj.a(2)+obj.p(2)+sqrt(sum((e_xy(2,:)-O_xy).^2)) obj.a(3)+obj.p(3)+sqrt(sum((e_xy(3,:)-O_xy).^2))];
            [k_l,l_l]=circcirc(obj.F(3,1),obj.F(3,2),obj.r(3),obj.F(2,1),obj.F(2,2),obj.r(2));
            [k_r,l_r]=circcirc(obj.F(3,1),obj.F(3,2),obj.r(3),obj.F(1,1),obj.F(1,2),obj.r(1));
            [k_t,l_t]=circcirc(obj.F(2,1),obj.F(2,2),obj.r(2),obj.F(1,1),obj.F(1,2),obj.r(1));
            [~,I1]=min(k_l);
            [~,I2]=max(k_r);
            [~,I3]=max(l_t);
            theta1=linspace(atan2d(l_r(I2)-obj.F(1,2),k_r(I2)-obj.F(1,1)),atan2d(l_t(I3)-obj.F(1,2),k_t(I3)-obj.F(1,1)),100);
            theta2=linspace(atan2d(l_t(I3)-obj.F(2,2),k_t(I3)-obj.F(2,1)),360+atan2d(l_l(I1)-obj.F(2,2),k_l(I1)-obj.F(2,1)),100);
            theta3=linspace(atan2d(l_l(I1)-obj.F(3,2),k_l(I1)-obj.F(3,1)),atan2d(l_r(I2)-obj.F(3,2),k_r(I2)-obj.F(3,1)),100);
            obj.ws_x=[obj.r(1)*cosd(theta1)+obj.F(1,1), obj.r(2)*cosd(theta2)+obj.F(2,1), obj.r(3)*cosd(theta3)+obj.F(3,1)];
            obj.ws_y=[obj.r(1)*sind(theta1)+obj.F(1,2), obj.r(2)*sind(theta2)+obj.F(2,2), obj.r(3)*sind(theta3)+obj.F(3,2)];
        end

%         function trace(obj,x_values,y_values)
%             for i=1:5
% 
%             end
%         end
    end
end
