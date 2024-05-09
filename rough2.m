clc
clear all
close all
warning off
format bank
%Step 1--Import Dataset
% Data=readtable('Mall_Customers.csv');
% %Step 2--Take only the AnnualIncome & Spedingscore columns(4 & 5)
% Data=Data(:,4:5);
% %Standardization
% Data.AnnualIncome=(Data.AnnualIncome-mean(Data.AnnualIncome))/std(Data.AnnualIncome);
% Data.SpendingScore=(Data.SpendingScore-mean(Data.SpendingScore))/std(Data.SpendingScore);
% %Convert the data from table to array 
% Data=table2array(Data);

Data = randn(800, 2);
x=Data(:,1);
y=Data(:,2);
scatter(x,y,'filled');
k=5;
cx=[-1 0 -1 3 -1.5];%initial mean of x
cy=[1 2 -1 2 1];%initial mean of y
mean_oldx=cx;
mean_newx=cx;
mean_oldy=cy;
mean_newy=cy;
outputx=cell(k,1);
outputy=cell(k,1);
temp=0;
    hold on;
    pause(3);
    scatter(cx,cy,'rx','LineWidth',20);
    hold off;
    pause(3);
while(temp==0)
    clf
    scatter(x,y,'filled');
    hold on;
    plot(cx,cy,'rx','LineWidth',20);
    mean_oldx=mean_newx;
    mean_oldy=mean_newy;
    for ij=1:length(x)
        mina=[];
        mu=x(ij);
        nu=y(ij);
     for mk=1:length(cx)
         mina=[mina sqrt((mu-cx(mk))^2+(nu-cy(mk))^2)];
     end
     [gc index]=min(mina);
     plot([x(ij),cx(index)],[y(ij),cy(index)],'Color',[0.75 0.75 0.75] );
         plot(cx,cy,'rx','LineWidth',20);
     outputx{index}=[outputx{index} mu];
     outputy{index}=[outputy{index} nu];
    end
    hold off;
    pause(0.5);
    gmckx=[];
    gmcky=[];
    for i=1:k
        gmckx=[gmckx mean(outputx{i})];
        gmcky=[gmcky mean(outputy{i})];
    end
    cx=gmckx;
    cy=gmcky;
    mean_newx=cx;
        mean_newy=cy;
        gum=0;
        bum=0;
    if(mean_newx==mean_oldx)
        gum=1;
    end
    if(mean_newy==mean_oldy)
        bum=1;
    end
    if(gum==1 && bum==1)
        temp=1;
    else
            outputx=cell(k,1);
            outputy=cell(k,1);
    end
end