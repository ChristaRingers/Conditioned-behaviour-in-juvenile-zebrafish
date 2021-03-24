classdef FunObj
    methods (Static)
        
        %%this are all the functions you need to run the above code
        function [ ITI,ITImean_std]=mean_duration_entrance(X,D,info,flipping,V)
            for i=1:size(X,1)
                for j=1:size(X,2)
                    x=[];x=X{i,j};d=[];d=cumsum(D{i,j});
                    v=[];v=V{i,j};
                    if (strcmp('NO',flipping))
                        border=info{3,j}
                        x(x<=border)=1;
                        x(x>border)=0;
                        xdiff=[];xdiff=diff(x);
                        ENTRANCE=[];ENTRANCE=find(xdiff==1);
                        EXIT=[];EXIT=find(xdiff==-1);
                        if (max(size(EXIT))>1)
                            if(d(EXIT(1))-d(ENTRANCE(1))<0)
                                EXIT(1)=[];
                            end
                            if(d(EXIT(end))-d(ENTRANCE(end))<0)
                                ENTRANCE(end)=[];
                            end
                            vel_tot=[];vel_mean=[];
                            for t=1:size(ENTRANCE)
                                vel_tot=[vel_tot,v(ENTRANCE(t):EXIT(t))];
                                vel_mean(t,1)=nanmean(v(ENTRANCE(t):EXIT(t)));
                            end
                            ITImean_std(i,j,4)=nanmean(vel_tot(:));
                            ITI(i,j)={[d(EXIT)'-d(ENTRANCE)']./1000};
                            ITImean_std(i,j,1)=nanmean(d(EXIT)'-d(ENTRANCE)')./1000;
                            ITImean_std(i,j,2)=nanmean(((d(EXIT)'-d(ENTRANCE)')./1000).*vel_mean);
                            ITImean_std(i,j,3)=size([d(EXIT)'-d(ENTRANCE)'],1);
                        else
                            ITI(i,j)={0};
                            ITImean_std(i,j,1)=0;
                            ITImean_std(i,j,2)=0;
                            ITImean_std(i,j,3)=0;
                            ITImean_std(i,j,4)=0;
                        end
                    else
                        if(i<flipping)
                            border=info{3,j}
                            x(x<=border)=1;
                            x(x>border)=0;
                        else
                            border=info{4,j}
                            x(x<border)=0;
                            x(x>=border)=1;
                        end
                        xdiff=[];xdiff=diff(x);
                        ENTRANCE=[];ENTRANCE=find(xdiff==1);
                        EXIT=[];EXIT=find(xdiff==-1);
                        if (max(size(EXIT))>1)
                            if(d(EXIT(1))-d(ENTRANCE(1))<0)
                                EXIT(1)=[];
                            end
                            if(d(EXIT(end))-d(ENTRANCE(end))<0)
                                ENTRANCE(end)=[];
                            end
                            vel_tot=[];vel_mean=[];
                            for t=1:size(ENTRANCE)
                                vel_tot=[A,v(ENTRANCE(t):EXIT(t))];
                                vel_mean(t,1)=nanmean(v(ENTRANCE(t):EXIT(t)));
                            end
                            ITImean_std(i,j,4)=nanmean(vel_tot(:));
                            ITI(i,j)={[d(EXIT)'-d(ENTRANCE)']./1000};
                            ITImean_std(i,j,1)=nanmean(d(EXIT)'-d(ENTRANCE)')./1000;
                            ITImean_std(i,j,2)=nanmean(((d(EXIT)'-d(ENTRANCE)')./1000).*vel_mean);
                            ITImean_std(i,j,3)=size([d(EXIT)'-d(ENTRANCE)'],1);
                        else
                            ITI(i,j)={0};
                            ITImean_std(i,j,1)=0;
                            ITImean_std(i,j,2)=0;
                            ITImean_std(i,j,3)=0;
                            ITImean_std(i,j,4)=0;
                        end
                    end
                end
            end
        end
        function [index]=LI_with_inversion_of_the_rule_baseline_refered(X_vector,Y_vector,D_vector,info_vector,change_rule)
            block=size(X_vector,1);
            fish_n=size(X_vector,2);
            learning_index=zeros(block,fish_n);
            for block_n=1:block
                for i=1: fish_n
                    
                    if(size(info_vector,1)==4)
                        if(block_n <change_rule)
                            threshold=info_vector{3,i};
                            threshold=cell2mat(info_vector(3,i));
                        else
                            if(isempty(info_vector{4,i}))
                                threshold=info_vector{3,i};
                                threshold=cell2mat(info_vector(3,i));
                            else
                                threshold=info_vector{4,i};
                                threshold=cell2mat(info_vector(4,i));
                            end
                        end
                    else
                        threshold=info_vector{3,i};
                        threshold=cell2mat(info_vector(3,i));
                    end
                    X=X_vector{block_n,i};
                    Y=Y_vector{block_n,i};
                    D=D_vector{block_n,i};
                    [timespenttemp,timetotal]=timespent(X,D,threshold);
                    if(strcmp(change_rule,'NO'))
                        if(block_n==1 )
                            timebasetemp(i)=((timespenttemp/timetotal));
                            index(block_n,i)=((timespenttemp/timetotal));
                        else
                            index(block_n,i)=(timebasetemp(i)-((timespenttemp/timetotal)))/timebasetemp(i);
                        end
                        
                    else
                        if((block_n==1)) %
                            timebasetemp(i)=((timespenttemp/timetotal));
                            index(block_n,i)=((timespenttemp/timetotal));
                        end
                        if(block_n>1)
                            if(block_n<change_rule)
                                index(block_n,i)=(timebasetemp(i)-((timespenttemp/timetotal)))/timebasetemp(i);
                            else
                                index(block_n,i)=((1-timebasetemp(i))-(((timetotal-(timespenttemp))/timetotal)))/(1-timebasetemp(i));
                                
                            end
                        end
                    end
                end
            end
            index(index<-1)=-1;
            clearvars X D Y threshold
        end
        function [Tspent, Ttotal]=timespent(X,D,threshold)
            timespent2=0;
            timetotal2=0;
            for j=2:length(X)
                if(X(j)<threshold)
                    timespent2 = timespent2 + D(j);
                end
                timetotal2= timetotal2 + D(j);
            end
            Tspent=timespent2;
            Ttotal=timetotal2;
        end
        function [xnew,ynew,dnew,pdf,vnew]=get_time_bins(X,Y,D,V,time_bin,size_arena,min_size_bin)%%timebin in min
            xnew=[];ynew=[];dnew=[];pdf=[];vnew=[];
            for i=1:size(X,2)
                size_arena2=max(size_arena{1:2,i});
                xnew2=cell([0]);ynew2=cell([0]);dnew2=cell([0]);pdf2=cell([0]);vnew2=cell([0]);
                index2=1;
                for j=1:(size(X,1))
                    x=X{j,i};
                    y=Y{j,i};
                    v=V{j,i};
                    if(size(x,2)==1)
                        v=v';
                    end
                    d=D{j,i};
                    if(size(x,2)==1)
                        x=x';
                    end
                    if(size(y,2)==1)
                        y=y';
                    end
                    if(size(d,2)==1)
                        d=d';
                    end
                    dtemp=[];
                    dtemp=diff(floor(cumsum(d)./(time_bin*60000)));
                    dtemp2=[];dtemp2=[1,find(dtemp)];
                    for index=1:size(dtemp2,2)
                        if(index<size(dtemp2,2))
                            xnew2(index2)={x(dtemp2(index):dtemp2(index+1))};
                            ynew2(index2)={y(dtemp2(index):dtemp2(index+1))};
                            vnew2(index2)={v(dtemp2(index):dtemp2(index+1))};
                            dnew2(index2)={d(dtemp2(index):dtemp2(index+1))};
                            pdf2(index2)={get_PDF(x(dtemp2(index):dtemp2(index+1))',y(dtemp2(index):dtemp2(index+1))',size_arena2)};
                            index2=index2+1;
                        else
                            if(max(size(x(dtemp2(index):end)))>(min_size_bin))
                                xnew2(index2)={x(dtemp2(index):end)};
                                ynew2(index2)={y(dtemp2(index):end)};
                                dnew2(index2)={d(dtemp2(index):end)};
                                vnew2(index2)={v(dtemp2(index):end)};
                                pdf2(index2)={get_PDF(x(dtemp2(index):end)',y(dtemp2(index):end)',size_arena2)};
                                index2=index2+1;
                            end
                        end
                    end
                    clearvars x d y v
                end
                if(size(xnew,1)>size(xnew2',1) && size(xnew,1)>0 )
                    xnew=[xnew,[xnew2';NaN(size(xnew,1)-size(xnew2',1))]];
                    ynew=[ynew,[ynew2';NaN(size(xnew,1)-size(xnew2',1))]];
                    dnew=[dnew,[dnew2';NaN(size(xnew,1)-size(xnew2',1))]];
                    pdf=[pdf,[pdf2';NaN(size(xnew,1)-size(xnew2',1))]];
                    vnew=[vnew,[vnew2';NaN(size(xnew,1)-size(xnew2',1))]];
                else
                    xnew=[xnew,xnew2'];ynew=[ynew,ynew2'];dnew=[dnew,dnew2'];pdf=[pdf,pdf2'];vnew=[vnew,vnew2'];
                end
            end
        end
        function varargout=shadedErrorBar(x,y,errBar,varargin)
            % generatcontinuous error bar area around a line plot
            %
            % function H=shadedErrorBar(x,y,errBar, ...)
            %
            % Purpose
            % Makes a 2-d line plot with a pretty shaded error bar made
            % using patch. Error bar color is chosen automatically.
            %
            %
            % Inputs (required)
            % x - vector of x values [optional, can be left empty]
            % y - vector of y values or a matrix of n observations by m cases
            %     where m has length(x);
            % errBar - if a vector we draw symmetric errorbars. If it has a size
            %          of [2,length(x)] then we draw asymmetric error bars with
            %          row 1 being the upper bar and row 2 being the lower bar
            %          (with respect to y). ** alternatively ** errBar can be a
            %          cellArray of two function handles. The first defines which
            %          statistic the line should be and the second defines the
            %          error bar.
            %
            % Inputs (optional, param/value pairs)
            % 'lineProps' - ['-k' by default] defines the properties of
            %             the data line. e.g.:
            %             'or-', or {'-or','markerfacecolor',[1,0.2,0.2]}
            % 'transparent' - [true  by default] if true, the shaded error
            %               bar is made transparent. However, for a transparent
            %               vector image you will need to save as PDF, not EPS,
            %               and set the figure renderer to "painters". An EPS
            %               will only be transparent if you set the renderer
            %               to OpenGL, however this makes a raster image.
            % 'patchSaturation'- [0.2 by default] The saturation of the patch color.
            %
            %
            %
            % Outputs
            % H - a structure of handles to the generated plot objects.
            %
            %
            % Examples:
            % y=randn(30,80);
            % x=1:size(y,2);
            %
            % 1)
            % shadedErrorBar(x,mean(y,1),std(y),'lineprops','g');
            %
            % 2)
            % shadedErrorBar(x,y,{@median,@std},'lineprops',{'r-o','markerfacecolor','r'});
            %
            % 3)
            % shadedErrorBar([],y,{@median,@(x) std(x)*1.96},'lineprops',{'r-o','markerfacecolor','k'});
            %
            % 4)
            % Overlay two transparent lines:
            % clf
            % y=randn(30,80)*10;
            % x=(1:size(y,2))-40;
            % shadedErrorBar(x,y,{@mean,@std},'lineprops','-r','transparent',1);
            % hold on
            % y=ones(30,1)*x; y=y+0.06*y.^2+randn(size(y))*10;
            % shadedErrorBar(x,y,{@mean,@std},'lineprops','-b','transparent',1);
            % hold off
            %
            %
            % Rob Campbell - November 2009
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Parse input arguments
            narginchk(3,inf)
            
            params = inputParser;
            params.CaseSensitive = false;
            params.addParameter('lineProps', '-k', @(x) ischar(x) | iscell(x));
            params.addParameter('transparent', true, @(x) islogical(x) || x==0 || x==1);
            params.addParameter('patchSaturation', 0.2, @(x) isnumeric(x) && x>=0 && x<=1);
            
            params.parse(varargin{:});
            
            %Extract values from the inputParser
            lineProps =  params.Results.lineProps;
            transparent =  params.Results.transparent;
            patchSaturation = params.Results.patchSaturation;
            
            if ~iscell(lineProps), lineProps={lineProps}; end
            
            
            %Process y using function handles if needed to make the error bar dynamically
            if iscell(errBar)
                fun1=errBar{1};
                fun2=errBar{2};
                errBar=fun2(y);
                y=fun1(y);
            else
                y=y(:).';
            end
            
            if isempty(x)
                x=1:length(y);
            else
                x=x(:).';
            end
            
            
            %Make upper and lower error bars if only one was specified
            if length(errBar)==length(errBar(:))
                errBar=repmat(errBar(:)',2,1);
            else
                s=size(errBar);
                f=find(s==2);
                if isempty(f), error('errBar has the wrong size'), end
                if f==2, errBar=errBar'; end
            end
            
            if length(x) ~= length(errBar)
                error('length(x) must equal length(errBar)')
            end
            
            
            %Log the hold status so we don't change
            initialHoldStatus=ishold;
            if ~initialHoldStatus, hold on,  end
            
            H = makePlot(x,y,errBar,lineProps,transparent,patchSaturation);
            
            if ~initialHoldStatus, hold off, end
            
            if nargout==1
                varargout{1}=H;
            end
            
            
            
            function H = makePlot(x,y,errBar,lineProps,transparent,patchSaturation)
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Plot to get the parameters of the line
                H.mainLine=plot(x,y,lineProps{:},'linewidth',2);
                
                
                % Work out the color of the shaded region and associated lines.
                % Here we have the option of choosing alpha or a de-saturated
                % solid colour for the patch surface.
                mainLineColor=get(H.mainLine,'color');
                edgeColor=mainLineColor+(1-mainLineColor)*0.55;
                
                if transparent
                    faceAlpha=patchSaturation;
                    patchColor=mainLineColor;
                else
                    faceAlpha=1;
                    patchColor=mainLineColor+(1-mainLineColor)*(1-patchSaturation);
                end
                
                
                %Calculate the error bars
                uE=y-errBar(1,:);
                lE=y+errBar(2,:);
                
                
                %Add the patch error bar
                
                
                
                %Make the patch
                yP=[lE,fliplr(uE)];
                xP=[x,fliplr(x)];
                
                %remove nans otherwise patch won't work
                xP(isnan(yP))=[];
                yP(isnan(yP))=[];
                
                
                H.patch=patch(xP,yP,1,'facecolor',patchColor, ...
                    'edgecolor','none', ...
                    'facealpha',faceAlpha);
                
                
                %Make pretty edges around the patch.
                H.edge(1)=plot(x,lE,'-','color',edgeColor);
                H.edge(2)=plot(x,uE,'-','color',edgeColor);
                
                
                
                uistack(H.mainLine,'top') % Bring the main line to the top
            end
            
            
        end
        function [pdf]=get_PDF(X,Y,dim_max)
            x= 0 : dim_max/25 : dim_max; % axis x, which you want to see
            y = 0 :dim_max/25: dim_max; % axis y, which you want to see
            [Xtemp,Ytemp] = meshgrid(x,y); % important for "surf" - makes defined grid
            pdf = hist3([X , Y],{x y}); % standard hist3 (calculated for yours axis
            pdf_normalize = (((pdf./ length(X)))); % normalization means devide it by length of
            pdf_normalize =(pdf_normalize)';% data_x (or data_y)
            pdf=pdf_normalize;
        end
        function [PDF_MEAN] = plot_heat_maps_averaged(PDF)
            x=size(PDF{1,1},1);
            y=size(PDF{1,1},2);
            for i=  1:size(PDF,1)
                for j = 1:size(PDF,2)
                    pdf(1:x,1:y,j)=PDF{i,j}
                end
                PDF_MEAN(i,1)={mean(pdf(1:end,1:end,:),3)};
            end
        end
        function [Binned_time]=new_binned_plot_with_flipping_color(X,Y,D,info_vector,color,time_bin,base_time,Flip_time)
            for i=1:size(X,2)
                max_size=max(info_vector{3,i},info_vector{2,i});
                for j=1:(size(X,1))
                    x=X{j,i};
                    y=Y{j,i};
                    d=D{j,i};
                    %         [x_3,y_3]=remove_pos_spike(x, y); you can implement here a
                    %         function to filter your data if a lot of detection errors are
                    %         present
                    %         x=[];  y=[];
                    %         x=x_3;
                    %         y=y_3;
                    if(size(x,2)==1)
                        x=x';
                    end
                    if(size(y,2)==1)
                        y=y';
                    end
                    if(size(d,2)==1)
                        d=d';
                    end
                    if(size(info_vector,1)==4)
                        if( strcmp(Flip_time,'NO'))
                            thresh=info_vector(3,:);
                            border=thresh{i};
                        else
                            if(j<Flip_time  )
                                thresh=info_vector(3,:);
                                border=thresh{i};
                            else
                                if(isempty(info_vector{4,1}))
                                    thresh=info_vector(3,:);
                                    border=thresh{i};
                                else
                                    thresh=info_vector(4,:);
                                    border=thresh{i};
                                end
                            end
                        end
                    else
                        thresh=info_vector(3,:);
                        border=thresh{i};
                    end
                    if(j==1)
                        time=0;
                        for k=max(size(x)):-1:1
                            
                            if(time<base_time && k>1)
                                time=time+d(k);
                            else
                                Binned_time(j,i)={divide_time_in_red(x(k:end),d(k:end),border,time_bin)};
                                break
                            end
                        end
                        clearvars k
                    else
                        if (strcmp(Flip_time,'NO'))
                            Binned_time(j,i)={divide_time_in_red(x,d,border,time_bin)};
                        else
                            if(j<Flip_time )
                                Binned_time(j,i)={divide_time_in_red(x,d,border,time_bin)};
                            else
                                Binned_time(j,i)={divide_time_in_gray(x,d,border,time_bin)};
                            end
                        end
                        clearvars x time  d y
                    end
                end
                clearvars time timetot point percentage_time
            end
        end
        function [] = plot_boundary(M, t)
            for s = 1:size(M,2)-1
                plot([M(s),M(s)], t, 'k','linewidth', 3);
                hold on
            end
            hold on
        end
        function [Mean_std]=plot_binned(Bin,color,binning)
            Xtot=NaN(size(Bin,2),size(Bin,1)*binning);
            for j=1:size(Bin,2)
                X=[];
                for i=1:size(Bin,1)
                    x=NaN(1,binning);
                    x=Bin{i,j};
                    X=[X,x];
                end
                plot(X)
                hold on
                size(X)
                X_tot(j,1:size(X,2))=X;
            end
            errorbar(1:size(X_tot,2),nanmean(X_tot,1),nanstd(X_tot,1)./sqrt(size(X_tot,1)),color ,'linewidth', 2 )
            xlim([0,size(X_tot,2)+1]);
            hold on
            Mean_std(1,1:size(X_tot,2))=nanmean(X_tot,1);
            Mean_std(2,1:size(X_tot,2))=nanstd(X_tot,1)./sqrt(size(X_tot,1));
        end
        function distance=calculate_distance(xpos,ypos)
            distance=(((xpos(2:end)-xpos(1:end-1)).^2)+((ypos(2:end)-ypos(1:end-1)).^2)).^(0.5);
        end
        function [Dist_raw,Distance_sum,Velocity,binned_vel,D_red]=binned_plot_velocity(X,Y,D,info_vector,Normalized,Flip_time,time_bin, base_time,V_vector)
            for i=1:size(X,2)
                max_size=max(info_vector{3,i},info_vector{2,i});
                for j=1:(size(X,1))
                    x=X{j,i};        y=Y{j,i};        d=D{j,i};
                    %         [x_3,y_3]=remove_pos_spike(x, y);if you want to filter your data
                    %         because of tracking errors implement it here
                    %         x=[];  y=[];    x=x_3;y=y_3;
                    if(size(x,2)==1)
                        x=x';
                    end
                    if(size(y,2)==1)
                        y=y';
                    end
                    if(size(d,2)==1)
                        d=d';
                    end
                    if(size(info_vector,1)==4)
                        if( strcmp(Flip_time,'NO'))
                            thresh=info_vector(3,:);
                            border=thresh{i};
                        else
                            if(j<Flip_time  )
                                thresh=info_vector(3,:);
                                border=thresh{i};
                            else
                                if(isempty(info_vector{4,i}))
                                    thresh=info_vector(3,:);
                                    border=thresh{i};
                                else
                                    thresh=info_vector(4,:);
                                    border=thresh{i};
                                end
                            end
                        end
                    else
                        thresh=info_vector(3,:);
                        border=thresh{i};
                    end
                    if(j==1)
                        time=0;
                        for k=max(size(x)):-1:1
                            if(time<base_time && k>1)
                                time=time+d(k);
                            else
                                break
                            end
                        end
                        Distance(j,i)=calculate_distance(x(k:end),y(k:end));
                        %             {(((x(k+1:end)-x(k:end-1)).^2)+((y(k+1:end)-y(k:end-1)).^2)).^(0.5)};
                        Dist_raw(j,i)={Distance{j,i}.*(120/max_size)};
                        Velocity=V_vector;
                        if (Normalized)
                            scale_factor=mean(Velocity{1,i});
                        else
                            scale_factor=1;
                        end
                        Velocity(j,i)={Velocity{j,i}./scale_factor};
                        v_temporary=[];v_temporary=Velocity{j,i};
                        binned_vel(j,i)={divide_velocity_bin(v_temporary(k+1:end),d(k+1:end),time_bin)};
                        binned_vel(j,i)={binned_vel{j,i}};
                        temp_dist_calc=[];
                        temp_dist_calc= Distance{j,i};
                        D_red(j,i)=sum(temp_dist_calc(x(k+1:end)<border)).*(0.12/max_size);
                        Distance_sum(j,i)=sum(Distance{j,i}).*(0.12/max_size);
                        clearvars k
                    else
                        Distance(j,i)=calculate_distance(x,y);
                        %             {(((x(2:end)-x(1:end-1)).^2)+((y(2:end)-y(1:end-1)).^2)).^(0.5)};
                        Dist_raw(j,i)={Distance{j,i}.*(120/max_size)};
                        Velocity(j,i)={Velocity{j,i}./scale_factor};
                        v_temporary=[];v_temporary=Velocity{j,i};
                        binned_vel(j,i)={divide_velocity_bin(v_temporary,d(2:end),time_bin)};
                        binned_vel(j,i)={binned_vel{j,i}};
                        if( strcmp(Flip_time,'NO'))
                            temp_dist_calc=[];
                            temp_dist_calc= Distance{j,i};
                            D_red(j,i)=sum(temp_dist_calc(x(1:end-1)<border)).*(0.12/max_size);
                        else
                            if(j<Flip_time  )
                                temp_dist_calc=[];
                                temp_dist_calc= Distance{j,i};
                                D_red(j,i)=sum(temp_dist_calc(x(1:end-1)<border)).*(0.12/max_size);
                            else
                                temp_dist_calc=[];
                                temp_dist_calc= Distance{j,i};
                                D_red(j,i)=sum(temp_dist_calc(x(1:end-1)>border)).*(0.12/max_size);
                            end
                        end
                        Distance_sum(j,i)=sum(Distance{j,i}).*(0.12/max_size);
                        clearvars x time  d y
                    end
                end
                clearvars time timetot point percentage_time
            end
            
        end
        function [velocity_bin]=divide_velocity_bin(V,D,time_bin)
            time_tot=D(1);
            point=1;
            counting=1;
            while(time_tot>=time_bin)
                velocity_bin(point)= 0;
                time_tot=abs(time_tot-time_bin);
                point=point+1;
            end
            for n=1:(length (V))
                time_tot=time_tot+D(n);
                if(time_tot<time_bin)
                    if(n==length(V) && time_tot > time_bin/20 )
                        velocity_bin(point)= mean(V(counting:n));
                        time_tot=0;
                        point=point+1;
                        counting=n;
                    end
                else
                    velocity_bin(point)= mean(V(counting:n));
                    time_tot/60000
                    time_tot=0;
                    point=point+1;
                    counting=n;
                end
            end
        end
        function [freezing_percentage,freezing_time]=freezing_new(X,Y,D,parameters,info_vector) % dist to freeze in mm , sec to freeze in sec
            dist_to_freeze=parameters(2);
            second_to_freeze=parameters(1);
            for j=1:size(X,1)
                for i=1:size(X,2)
                    tempd=[];tempd= D{j,i};
                    tempx= [];
                    tempy= [];
                    tempx= X{j,i};
                    tempy= Y{j,i};
                    deltax=[];
                    deltay=[];
                    dist_tot=[];
                    deltax=diff(tempx);
                    deltay=diff(tempy);
                    dist_tot=sqrt((deltax.^2 )+(deltay.^2 ))./(120/min(info_vector{1:2,i}));
                    dist_cumulative=floor(cumsum(dist_tot)./dist_to_freeze);
                    dist_cum=diff(dist_cumulative);
                    moved=[];moved=[find(dist_cum>0)];
                    deltat=[];
                    deltat(1)=sum(tempd(1:moved(1))) ;
                    for t= 1: max(size(moved))
                        if(t< max(size(moved)))
                            deltat(t)=sum(tempd((moved(t)+1):moved(t+1))) ;
                        else
                            deltat(t)=sum(tempd((moved(t)+1):moved(end))) ;
                        end
                    end
                    freezing_time(j,i)=sum(deltat(deltat>=(second_to_freeze*1000)));
                    freezing_percentage(j,i)=freezing_time(j,i)./sum(deltat);
                end
            end
            
        end
        function [percentage_time]=divide_time_in_red(X,D,border,time_bin) %%calculate the time on the left side of the arena
            time_tot=0;
            time=0;
            point=1;
            sum(D)/60000
            for n=1:(length (X))
                if(time_tot<=time_bin)
                    time_tot=time_tot+D(n);
                    if(X(n)<border)
                        time=time+D(n);
                    end
                    if(n==length(X) && time_tot > 60000 )
                        percentage_time(point)= time/time_tot;
                        time=0;
                        time_tot=0;
                        point=point+1;
                        continue
                    else
                        if(n==length(X) && time_tot <= 60000 && time_tot>30000)
                            percentage_time(point)=percentage_time(point-1);
                            time=0;
                            time_tot=0;
                            point=point+1;
                            continue
                        end
                    end
                else
                    percentage_time(point)= time/time_tot;
                    time_tot/60000
                    time=0;
                    time_tot=0;
                    point=point+1;
                end
            end
        end
        function [binned_dist_red,binned_dist_all] = dist_in_red_time_bin(X,Y,D,info_vector,Normalized,Flip_time,time_bin,base_time,Normalise_to_base,czcs)
            
            for i=1:size(X,2)
                
                max_size=max(info_vector{3,i},info_vector{2,i});
                for j=1:(size(X,1))
                    
                    x=X{j,i};
                    y=Y{j,i};
                    d=D{j,i};
                    %         [x_3,y_3]=remove_pos_spike(x, y); to filter your data
                    %         x=[];  y=[];
                    %
                    %         x=x_3;
                    %         y=y_3;
                    if(size(x,2)==1)
                        x=x';
                    end
                    if(size(y,2)==1)
                        y=y';
                    end
                    if(size(d,2)==1)
                        d=d';
                    end
                    if(size(info_vector,1)==4)
                        if( strcmp(Flip_time,'NO'))
                            thresh=info_vector(3,:);
                            border=thresh{i};
                        else
                            if(j<Flip_time  )
                                thresh=info_vector(3,:);
                                border=thresh{i};
                            else
                                if(isempty(info_vector{4,1}))
                                    thresh=info_vector(3,:);
                                    border=thresh{i};
                                else
                                    thresh=info_vector(4,:);
                                    border=thresh{i};
                                end
                            end
                        end
                    else
                        thresh=info_vector(3,:);
                        border=thresh{i};
                    end
                    if(j==1)
                        time=0;
                        for k=max(size(x)):-1:1
                            if(time<base_time && k>1)
                                time=time+d(k);
                            else
                                break
                            end
                        end
                        Distance(j,i)= {(((x(k+1:end)-x(k:end-1)).^2)+((y(k+1:end)-y(k:end-1)).^2)).^(0.5)};
                        Dist_raw(j,i)={Distance{j,i}.*(12/max_size)};
                        [binned_dist_red(j,i),binned_dist_all(j,i)]=divide_distance_bin( d(k+1:end),Dist_raw{j,i},x(k+1:end),time_bin,border,1,czcs);
                        if(Normalise_to_base)
                            normalise_factor(i)=mean(binned_dist_red{j,i});
                        else
                            normalise_factor(i)=1;
                        end
                        binned_dist_red(j,i)={ binned_dist_red{j,i}./normalise_factor(i)};
                        clearvars k
                    else
                        Distance(j,i)= {(((x(2:end)-x(1:end-1)).^2)+((y(2:end)-y(1:end-1)).^2)).^(0.5)};
                        Dist_raw(j,i)={Distance{j,i}.*(12/max_size)};
                        if( strcmp(Flip_time,'NO'))
                            [binned_dist_red(j,i),binned_dist_all(j,i)]=divide_distance_bin(d,Dist_raw{j,i},x(2:end),time_bin,border,1,czcs);
                            binned_dist_red(j,i)={ binned_dist_red{j,i}./normalise_factor(i)};
                        else
                            if(j<Flip_time  )
                                [binned_dist_red(j,i),binned_dist_all(j,i)]=divide_distance_bin(d,Dist_raw{j,i},x(2:end),time_bin,border,1,czcs);
                                TT=[];TT= binned_dist_red(j,i);
                                binned_dist_red(j,i)={ binned_dist_red{j,i}./normalise_factor(i)};
                                if(j==Flip_time-1 && Normalise_to_base)
                                    normalise_factor(i)=mean(TT{:}(end-5:end));
                                end
                            else
                                [binned_dist_red(j,i),binned_dist_all(j,i)]=divide_distance_bin(d,Dist_raw{j,i},x(2:end),time_bin,border,2,czcs);
                                binned_dist_red(j,i)={ binned_dist_red{j,i}./normalise_factor(i)};
                            end
                        end
                        clearvars x time  d y
                    end
                end
                clearvars time timetot point percentage_time
            end
        end
        function [Mean_std,X_tot]=plot_binned_dist(Bin,color,binning)
            Xtot=NaN(size(Bin,2),size(Bin,1)*binning);
            for j=1:size(Bin,2)
                X=[];
                for i=1:size(Bin,1)
                    x=NaN(1,binning);
                    x=Bin{i,j};
                    X=[X,x];
                end
                plot(X)
                hold on
                size(X)
                X_tot(j,1:size(X,2))=X;
            end
            errorbar(1:size(X_tot,2),nanmean(X_tot,1),nanstd(X_tot,1)./sqrt(size(X_tot,1)),color ,'linewidth', 2 )
            xlim([0,size(X_tot,2)+1]);
            hold on
            Mean_std(1,1:size(X_tot,2))=nanmean(X_tot,1);
            Mean_std(2,1:size(X_tot,2))=nanstd(X_tot,1)./sqrt(size(X_tot,1));
        end
        function [dist_R,dist_A]=divide_distance_bin(time,D,X,time_bin,threshold,rule,czcs)% if rule is 1 x<threshold else x>threshold
            time_tot=0;
            point=1;
            counting=1;
            for n=1:(length (D))
                
                
                
                if(time_tot<time_bin)
                    time_tot=time_tot+time(n);
                    
                    if(n==length(D) && time_tot >(0.5*time_bin ))
                        if(rule==1)
                            tempD=[];tempD=D(counting:n);
                            dist_red(point)= sum(tempD(X(counting:n)<threshold));
                            dist_all(point)= sum(D(counting:n));
                            
                            time_tot=0;
                            point=point+1;
                            counting=n;
                        else
                            tempD=[];tempD=D(counting:n);
                            dist_red(point)= sum(tempD(X(counting:n)>threshold));
                            dist_all(point)= sum(D(counting:n));
                            
                            time_tot=0;
                            point=point+1;
                            counting=n;
                        end
                    end
                else
                    if(rule==1)
                        tempD=[];tempD=D(counting:n);
                        dist_red(point)= sum(tempD(X(counting:n)<threshold));
                        dist_all(point)= sum(D(counting:n));
                        time_tot/60000
                        time_tot=0;
                        point=point+1;
                        counting=n;
                    else
                        tempD=[];tempD=D(counting:n);
                        dist_red(point)= sum(tempD(X(counting:n)>threshold));
                        dist_all(point)= sum(D(counting:n));
                        time_tot/60000
                        time_tot=0;
                        point=point+1;
                        counting=n;
                    end
                end
            end
            dist_R= {dist_red};
            temp=[];
            if (czcs)
                temp=dist_red./(dist_all-dist_red);
                temp(temp>5)=5;
            else
                temp=(dist_all-dist_red);
            end
            dist_A={temp};
        end
        function [percentage_time]=get_mean_pos_in_timebin(X,D,border,time_bin,mean_or_median)
            time_tot=0;
            point=1;
            sum(D)/60000
            x_temporary=[];
            point_t=1;
            for n=1:(length (X))
                if(time_tot<time_bin)
                    time_tot=time_tot+D(n);
                    if(n==length(X) && time_tot > 60000 )
                        if(strcmp(mean_or_median,'mean'))
                            x_temporary= border-X(point_t:n);
                            percentage_time(point)= mean(x_temporary);
                        else
                            x_temporary= border-X(point_t:n);
                            percentage_time(point)= median(x_temporary);
                        end
                        x_temporary=[];
                        point_t=n;
                        time_tot=0;
                        point=point+1;
                    end
                else
                    if(strcmp(mean_or_median,'mean'))
                        x_temporary= border-X(point_t:n);
                        percentage_time(point)= mean(x_temporary);
                    else
                        x_temporary= border-X(point_t:n);
                        percentage_time(point)= median(x_temporary);
                    end
                    x_temporary=[];
                    time_tot=0;
                    point=point+1;
                    point_t=n;
                end
            end
        end
        function [Binned_pos]=binned_mean_position_flipping_color(X,Y,D,info_vector,color,time_bin,base_time,Flip_time,mean_or_median)
            for i=1:size(X,2)
                max_size=max(info_vector{3,i},info_vector{2,i});
                for j=1:(size(X,1))
                    x=X{j,i};
                    y=Y{j,i};
                    d=D{j,i};
                    %         [x_3,y_3]=remove_pos_spike(x, y);filter your data in case here
                    %         x=[];  y=[];
                    %         x=x_3;
                    %         y=y_3;
                    if(size(x,2)==1)
                        x=x';
                    end
                    if(size(y,2)==1)
                        y=y';
                    end
                    if(size(d,2)==1)
                        d=d';
                    end
                    if(size(info_vector,1)==4)
                        if( strcmp(Flip_time,'NO'))
                            thresh=info_vector(3,:);
                            border=thresh{i};
                        else
                            if(j<Flip_time  )
                                thresh=info_vector(3,:);
                                border=thresh{i};
                            else
                                if(isempty(info_vector{4,1}))
                                    thresh=info_vector(3,:);
                                    border=thresh{i};
                                else
                                    thresh=info_vector(4,:);
                                    border=thresh{i};
                                end
                            end
                        end
                    else
                        thresh=info_vector(3,:);
                        border=thresh{i};
                    end
                    if(j==1)
                        time=0;
                        for k=max(size(x)):-1:1
                            if(time<base_time && k>1)
                                time=time+d(k);
                            else
                                break
                            end
                        end
                        Binned_pos(j,i)={get_mean_pos_in_timebin(x(k:end),d(k:end),border,time_bin,mean_or_median)};
                        clearvars k
                    else
                        if (strcmp(Flip_time,'NO'))
                            Binned_pos(j,i)={get_mean_pos_in_timebin(x,d,border,time_bin,mean_or_median)};
                        else
                            if(j<Flip_time )
                                Binned_pos(j,i)={get_mean_pos_in_timebin(x,d,border,time_bin,mean_or_median)};
                            else
                                Binned_pos(j,i)={get_mean_pos_in_timebin(x,d,border,time_bin,mean_or_median)};
                            end
                        end
                        clearvars x time  d y
                    end
                end
                clearvars time timetot point percentage_time
            end
        end
        function [Mean_std]=plot_binned_vel(Bin,color)
            Bin
            X_tot=NaN(size(Bin,2),29);
            for j=1:size(Bin,2)
                X=[];
                for i=1:size(Bin,1)
                    x=Bin{i,j};
                    X=[X,x];
                end
                plot(X)
                hold on
                size(X)
                X_tot(j,1:size(X,2))=X;
            end
            errorbar(1:size(X_tot,2),nanmean(X_tot,1),nanstd(X_tot,1)./sqrt(size(X_tot,1)),color ,'linewidth', 2 )
            hold on
            Mean_std(1,1:size(X_tot,2))=nanmean(X_tot,1);
            Mean_std(2,1:size(X_tot,2))=nanstd(X_tot,1)./sqrt(size(X_tot,1));
        end
        function []=scattercloud(X,Y,filled,color,size)% 1 if you want the scatterplot filled
            data=-0.2:0.02:+0.2;
            for i=1:length(X)
                randomnum= datasample(data,1); % returns k observations sampled uniformly at random, with replacement, from the data in data.
                X(i)=X(i)+randomnum;
            end
            if(strcmp(filled,'filled'))
                scatter(X,Y,size,filled,color);
            else
                scatter(X,Y,size,color);
            end
        end
        
        function [ANGLE_ALL,V,V_raw,Peaks,MeanP,BASE_V,BASE_V_raw,PRE_Peak_response,Peak_response]=plot_response_to_only_shock_normalized(V_input,X,Y,info,window,switch_rule,internal_window,variability_threshold,division_in_time,time_block)
            % %%the varianle angles has the following information :
            % angles(line,1) angle of aproaching the midline
            % angles(line,2) x position of encounter of the midline
            % angles(line,3) x position "time_window" frames before meeting the midline
            % angles(line,4) y position of encounter of the midline
            % angles(line,5) y position "time_window" frames before meeting the midline
            % angles(line,6) mean velocity of appraching the midline
            
            % angles(line,7) angle of leaving the midline
            % angles(line,8) x position of encounter of the midline
            % angles(line,9) x position "time_window" frames after meeting the midline
            % angles(line,10) y position of encounter of the midline
            % angles(line,11) y position "time_window" frames after meeting the midline
            % angles(line,12) mean velocity of leaving the midline
            
            for i=1:size(X,2)
                for j=1:size(X,1)
                    line=1;
                    x=[];y=[];v=[];
                    x=X{j,i};
                    x=x(((end/division_in_time)*(time_block-1))+1:((end/division_in_time)*(time_block)));
                    v=V_input{j,i};
                    y=Y{j,i};
                    if( 1~=(strcmp(switch_rule,'NO')))
                        if(j>=switch_rule && size(info,2)==4 )
                            if(isempty(info{4,i}))
                                tresh=info{3,i}-variability_threshold;
                            else
                                tresh=info{4,i}-variability_threshold;
                            end
                        else
                            tresh=info{3,i}+variability_threshold;
                            
                        end
                    else
                        tresh=info{3,i}+variability_threshold;
                    end
                    if (j==1 && size(x,2)>28000)%% this is a control for the time division of the session which in our case were 30 min long
                        x=x(27000:end);
                        y=y(27000:end);
                        v=v(27000:end);
                    end
                    count=window;
                    A=[];
                    B=[];
                    angles=[];
                    for t=window:(max(size(x))-(window+1))
                        if(t>count && t<=(max(size(v))-(window+1)))
                            if( 1~=(strcmp(switch_rule,'NO')))
                                if(j<switch_rule )
                                    if (x(t)<tresh && x(t-internal_window)>tresh  && x(t+internal_window)<tresh)
                                        count=t+window-1;
                                        A(line,1:window)=v(t:t+window-1);
                                        if(t>window)
                                            HYP=[];CAT1=[];CAT2=[];
                                            CAT1=x(t+window)-x(t);
                                            CAT2=y(t+window)-y(t);
                                            HYP=sqrt(CAT1.^2+CAT2.^2);
                                            angles(line,8)=x(t);
                                            angles(line,9)=x(t+window);
                                            angles(line,10)=y(t);
                                            angles(line,11)=y(t+window);
                                            angles(line,12)=mean(mean(v(t:(t+window))));
                                            if(HYP==0)
                                                angles(line,7)=361;
                                                angles(line,12)=-0.1;
                                            else
                                                if(rad2deg(asin(CAT2/HYP))>=0)
                                                    angles(line,7)= ( rad2deg(acos(CAT1/HYP)));%post crossing thresh angle
                                                else
                                                    angles(line,7)= 360+(-( rad2deg(acos(CAT1/HYP))));%post crossing thresh angle
                                                end
                                            end
                                            HYP=[];CAT1=[];CAT2=[];
                                            CAT1=x(t)-x(t-window);
                                            CAT2=y(t)-y(t-window);
                                            HYP=sqrt(CAT1.^2+CAT2.^2);
                                            angles(line,2)=x(t);
                                            angles(line,3)=x(t-window);
                                            angles(line,4)=y(t);
                                            angles(line,5)=y(t-window);
                                            angles(line,6)=mean(mean(v((t-window):t)));
                                            if(HYP==0)
                                                angles(line,1)=361;
                                            else
                                                if(rad2deg(asin(CAT2/HYP))>=0)
                                                    angles(line,1)= ( rad2deg(acos(CAT1/HYP)));%pre crossing thresh angle
                                                else
                                                    angles(line,1)= 360+(-( rad2deg(acos(CAT1/HYP))));%pre crossing thresh angle
                                                end
                                            end
                                            B(line,1:window)=v(t-(window-1):t);
                                        end
                                        line=line+1;
                                    end
                                else
                                    if (x(t)>tresh && x(t-internal_window)<tresh && x(t+internal_window)>tresh)
                                        count=t+window-1;
                                        A(line,1:window)=v(t:t+window-1);
                                        if(t>window)
                                            HYP=[];CAT1=[];CAT2=[];
                                            CAT1=x(t+window)-x(t);
                                            CAT2=y(t+window)-y(t);
                                            HYP=sqrt(CAT1.^2+CAT2.^2);
                                            angles(line,8)=x(t);
                                            angles(line,9)=x(t+window);
                                            angles(line,10)=y(t);
                                            angles(line,11)=y(t+window);
                                            angles(line,12)=mean(mean(v(t:(t+window))));
                                            if(HYP==0)
                                                angles(line,7)=361;
                                                angles(line,12)=-0.1;
                                            else
                                                if(rad2deg(asin(CAT2/HYP))>=0)
                                                    angles(line,7)= ( rad2deg(acos(CAT1/HYP)));%post crossing thresh angle
                                                else
                                                    angles(line,7)= 360+(-( rad2deg(acos(CAT1/HYP))));%post crossing thresh angle
                                                end
                                            end
                                            HYP=[];CAT1=[];CAT2=[];
                                            CAT1=x(t)-x(t-window);
                                            CAT2=y(t)-y(t-window);
                                            HYP=sqrt(CAT1.^2+CAT2.^2);
                                            angles(line,2)=x(t);
                                            angles(line,3)=x(t-window);
                                            angles(line,4)=y(t);
                                            angles(line,5)=y(t-window);
                                            angles(line,6)=mean(mean(v((t-window):t)));
                                            if(HYP==0)
                                                angles(line,1)=361;
                                                angles(line,6)=-0.1;
                                            else
                                                if(rad2deg(asin(CAT2/HYP))>=0)
                                                    angles(line,1)= ( rad2deg(acos(CAT1/HYP)));%pre crossing thresh angle
                                                else
                                                    angles(line,1)= 360+(-( rad2deg(acos(CAT1/HYP))));%pre crossing thresh angle
                                                end
                                            end
                                            B(line,1:window)=v(t-(window-1):t);
                                        end
                                        line=line+1;
                                    end
                                end
                            else
                                if (x(t)<tresh && x(t-internal_window)>tresh && x(t+internal_window)<tresh)
                                    count=t+window-1;
                                    A(line,1:window)=v(t:t+window-1);
                                    if(t>window )
                                        HYP=[];CAT1=[];CAT2=[];
                                        CAT1=x(t+window)-x(t);
                                        CAT2=y(t+window)-y(t);
                                        HYP=sqrt(CAT1.^2+CAT2.^2);
                                        angles(line,8)=x(t);
                                        angles(line,9)=x(t+window);
                                        angles(line,10)=y(t);
                                        angles(line,11)=y(t+window);
                                        angles(line,12)=mean(mean(v(t:(t+window))));
                                        t+window
                                        if(HYP==0)
                                            angles(line,7)=361;
                                            angles(line,12)=-0.1;
                                        else
                                            if(rad2deg(asin(CAT2/HYP))>=0)
                                                angles(line,7)= ( rad2deg(acos(CAT1/HYP)));%post crossing thresh angle
                                            else
                                                angles(line,7)= 360+(-( rad2deg(acos(CAT1/HYP))));%post crossing thresh angle
                                            end
                                        end
                                        HYP=[];CAT1=[];CAT2=[];
                                        CAT1=x(t)-x(t-window);
                                        CAT2=y(t)-y(t-window);
                                        HYP=sqrt(CAT1.^2+CAT2.^2);
                                        angles(line,2)=x(t);
                                        angles(line,3)=x(t-window);
                                        angles(line,4)=y(t);
                                        angles(line,5)=y(t-window);
                                        angles(line,6)=mean(mean(v((t-window):t)));
                                        if(HYP==0)
                                            angles(line,1)=361;
                                            angles(line,6)=-0.1;
                                        else
                                            if(rad2deg(asin(CAT2/HYP))>=0)
                                                angles(line,1)= ( rad2deg(acos(CAT1/HYP)));%pre crossing thresh angle
                                            else
                                                angles(line,1)= 360+(-( rad2deg(acos(CAT1/HYP))));%pre crossing thresh angle
                                            end
                                        end
                                        B(line,1:window)=v(t-(window-1):t);
                                    end
                                    line=line+1;
                                end
                            end
                        end
                    end
                    ANGLE_ALL(j,i)={angles};
                    if(size(A,1)>1)
                        V(j,i)={mean(A,1)};
                        BASE_V(j,i)={mean(B,1)};
                        V_raw(j,i)={A};
                        BASE_V_raw(j,i)={B};
                        A=[];
                        B=[];
                        pks = mean( V{j,i});
                        Peaks(j,i)={pks};
                        pks_pre = mean( BASE_V{j,i});
                        Peaks_pre(j,i)={pks_pre};
                        Peak_response(j,i)=mean(pks);
                        PRE_Peak_response(j,i)=mean(pks_pre);
                        MeanP(j,i)=(mean(pks)./mean(pks_pre));
                        pks=[];locs=[];w=[];p=[];T=[];
                        pks_pre=[];locs_pre=[];w_pre=[];p_pre=[];T_pre=[];
                        hold on
                    else
                        if(size(A,1)==1)
                            V(j,i)={A};
                            BASE_V(j,i)={B};
                            V_raw(j,i)={A};
                            BASE_V_raw(j,i)={B};
                            A=[];
                            B=[];
                            pks = ( V{j,i});
                            Peaks(j,i)={pks};
                            pks_pre = ( BASE_V{j,i});
                            Peaks_pre(j,i)={pks_pre};
                            Peak_response(j,i)=mean(pks);
                            PRE_Peak_response(j,i)=mean(pks_pre);
                            MeanP(j,i)=(mean(pks)./mean(pks_pre));
                            pks=[];locs=[];w=[];p=[];T=[];
                            pks_pre=[];locs_pre=[];w_pre=[];p_pre=[];T_pre=[];
                            
                            hold on
                        else
                            V(j,i)={-0.1};
                            BASE_V(j,i)={-0.1};
                            V_raw(j,i)={A};
                            BASE_V_raw(j,i)={B};
                            A=[];
                            B=[];
                            pks = ( V{j,i});
                            Peaks(j,i)={pks};
                            pks_pre = ( BASE_V{j,i});
                            Peaks_pre(j,i)={pks_pre};
                            Peak_response(j,i)=(pks);
                            PRE_Peak_response(j,i)=(pks_pre);
                            MeanP(j,i)=((pks)./(pks_pre));
                            pks=[];locs=[];w=[];p=[];T=[];
                            pks_pre=[];locs_pre=[];w_pre=[];p_pre=[];T_pre=[];
                            hold on
                        end
                    end
                end
            end
        end
        function [mu ul ll] = circ_mean(alpha, w)
            %
            % mu = circ_mean(alpha, w)
            %   Computes the mean direction for circular data.
            %
            %   Input:
            %     alpha	sample of angles in radians
            %     [w		weightings in case of binned angle data]
            %
            %   Output:
            %     mu		mean direction
            %     ul    upper 95% confidence limit
            %     ll    lower 95% confidence limit
            %
            % PHB 7/6/2008
            %
            % References:
            %   Statistical analysis of circular data, N. I. Fisher
            %   Topics in circular statistics, S. R. Jammalamadaka et al.
            %   Biostatistical Analysis, J. H. Zar
            %
            % Circular Statistics Toolbox for Matlab
            
            % By Philipp Berens, 2009
            % berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
            
            % check vector size
            if size(alpha,2) > size(alpha,1)
                alpha = alpha';
            end
            
            if nargin<2
                % if no specific weighting has been specified
                % assume no binning has taken place
                w = ones(size(alpha));
            else
                if size(w,2) > size(w,1)
                    w = w';
                end
            end
            
            % compute weighted sum of cos and sin of angles
            r = w'*exp(1i*alpha);
            
            % obtain mean by
            mu = angle(r);
            
            % confidence limits if desired
            if nargout > 1
                t = circ_confmean(alpha,0.05,w);
                ul = mu + t;
                ll = mu - t;
            end
        end
        
    end
end
