function varargout = AFALA(varargin)
% AFALA M-file for AFALA.fig
%      AFALA, by itself, creates a new AFALA or raises the existing
%      singleton*.
%
%      H = AFALA returns the handle to a new AFALA or the handle to
%      the existing singleton*.
%
%      AFALA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AFALA.M with the given input arguments.
%
%      AFALA('Property','Value',...) creates a new AFALA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before calcLocalizationRandomGraph_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AFALA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help AFALA
% 2014-01-01 Flip Ambiguities with Transmission Range
% 2014-03-15 Error
% 2014-04-02 Polygon
% Last Modified by GUIDE v2.5 12-Apr-2018 23:43:27
% 颜色 r 红 g 绿 b 蓝 c 蓝绿 m 紫红 y 黄 k 黑 w 白

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AFALA_OpeningFcn, ...
                   'gui_OutputFcn',  @AFALA_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before AFALA is made visible.
function AFALA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AFALA (see VARARGIN)

% Choose default command line output for AFALA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AFALA wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AFALA_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- distance between vertex procedure.
function dSqr = distVertexes(x1, y1, x2, y2)
% x1, y1 is the cordinator of one vertex, x2, y2 is another.

deltax = x1 - x2;
deltay = y1 - y2;
dSqr = deltax.^2 + deltay.^2;
                        

%--judge whether triangle
function Triangle = IsTriangle(ab, ac, bc)
% ab, ac and bc is three edges^2

if ((sqrt(ab)+sqrt(ac)) > sqrt(bc)) && ((sqrt(ab)+sqrt(bc)) > sqrt(ac)) && ((sqrt(ac)+sqrt(bc)) > sqrt(ab)),
   Triangle = 1;
else
   Triangle = 0;
end


% --- calculate the interarea by criterion1-trilateration
function interarea = GetCircleIntersecByCriterion1(xo, yo, xp, yp, xq, yq, dop, doq, dpq, doi, dpi, dqi, rangeError, rangeVelocity)
disp(['xo:', num2str(xo),' yo:', num2str(yo),' xp:', num2str(xp),' yp:', num2str(yp),' xq:', num2str(xq),' yq:', num2str(yq),' dop:', num2str(dop),' doq:', num2str(doq),' dpq:', num2str(dpq),' doi:', num2str(doi),' dpi:', num2str(dpi),' dqi:', num2str(dqi),' rangeError:', num2str(rangeError)]);

%把速度和方向看成另一种误差因素参与运算
rangeError = rangeError + rangeVelocity;

if rangeError == 0,
    %the distance between so and sp
    r = sqrt((xp-xo)^2+(yp-yo)^2);%圆心距
    %如果两圆（oi与pi）相交
    if (doi+dpi>r) && abs(dpi-doi)<r,
        disp('相交1');
        %计算两圆的交点
        [xi,yi]=circcirc(xo, yo, doi, xp, yp, dpi);
        data = sortrows(double(vpa([xi;yi]', 5)),2)';
        xi1 = data(1,2); yi1 = data(2,2);
        xi2 = data(1,1); yi2 = data(2,1);
        %计算两个交点与sq的距离
        dqi1 = sqrt((xq-xi1)^2+(yq-yi1)^2);
        dqi2 = sqrt((xq-xi2)^2+(yq-yi2)^2);
        %确定哪一个交点是正确的
        if abs(dqi-dqi1) < 0.0001,
            interarea.area1 = 1;
            interarea.area2 = 0;
            return;
        elseif abs(dqi-dqi2) < 0.0001,
            interarea.area1 = 0;
            interarea.area2 = 1;
            return;
        else
            interarea.area1 = -1;
            interarea.area2 = -1;
            return;
        end
    else
        disp('两圆不相交1');
        interarea.area1 = -1;
        interarea.area2 = -1;
        return;
    end
else
    mesh_density=rangeError/3;
    
    %get the possible region boder of sq
    %大半圆环
    xoq1 = -(doq+rangeError):mesh_density:(doq+rangeError);%out
    yoq1 = sqrt((doq+rangeError).^2 - xoq1.^2);
    %小半圆环
    xoq2 = -(doq-rangeError):mesh_density:(doq-rangeError);%in
    yoq2 = sqrt((doq-rangeError).^2 - xoq2.^2);
    %线段1
    xoq31 = -(doq+rangeError):mesh_density:-(doq-rangeError);
    yoq31 = zeros(1,length(xoq31));
    %线段2
    xoq32 = (doq-rangeError):mesh_density:(doq+rangeError);
    yoq32 = zeros(1,length(xoq32));
    %边缘线
    Xoq1 = [xoq1, xoq2, xoq31, xoq32];
    Yoq1 = [yoq1, yoq2, yoq31, yoq32];
%     plot(Xoq1, Yoq1, 'b*');

    %左大半圆环
    xpq11 = (dop-rangeError)-(dpq+rangeError):mesh_density:dop;
    ypq11 = sqrt((dpq+rangeError).^2 - (xpq11-(dop-rangeError)).^2);
    %右大半圆环
    xpq12 = dop:mesh_density:(dop+rangeError)+(dpq+rangeError);
    ypq12 = sqrt((dpq+rangeError).^2 - (xpq12-((dop+rangeError))).^2);
    %左小半圆环
    xpq21 = (dop+rangeError)-(dpq-rangeError):mesh_density:dop;
    ypq21 = sqrt((dpq-rangeError).^2 - (xpq21-(dop+rangeError)).^2);
    %右小半圆环
    xpq22 = dop:mesh_density:(dop-rangeError)+(dpq-rangeError);
    ypq22 = sqrt((dpq-rangeError).^2 - (xpq22-(dop-rangeError)).^2);
    %线段1
    xpq31 = (dop-rangeError)-(dpq+rangeError):mesh_density:(dop+rangeError)-(dpq-rangeError);
    ypq31 = zeros(1,length(xpq31));
    %线段2
    xpq32 = (dop-rangeError)+(dpq-rangeError):mesh_density:(dop+rangeError)+(dpq+rangeError);
    ypq32 = zeros(1,length(xpq32));
    %线段2
    Xpq1 = [xpq11, xpq12, xpq21, xpq22, xpq31, xpq32];
    Ypq1 = [ypq11, ypq12, ypq21, ypq22, ypq31, ypq32];
%     plot(Xpq1, Ypq1, 'b*');
    
    %Define a circle xo, yo, doq-rangeError，外
    t = 0:0.01:2*pi;
    a = (doq-rangeError)*sin(t)+xo;
    b = (doq-rangeError)*cos(t)+yo;
    in_circle12_1 = inpolygon(Xpq1,Ypq1,a,b);
    index12_1 = find(in_circle12_1==0);
    
    %Define a circle xo, yo, doq+rangeError，内
    t = 0:0.01:2*pi;
    a = (doq+rangeError)*sin(t)+xo;
    b = (doq+rangeError)*cos(t)+yo;
    in_circle11_1 = inpolygon(Xpq1,Ypq1,a,b);
    index11_1 = find(in_circle11_1>0);

    %边缘交点
    Xpq1 = Xpq1(intersect(index12_1, index11_1));%a part of region border
    Ypq1 = Ypq1(intersect(index12_1, index11_1));
%     plot(Xpq1, Ypq1, 'b*');

    %Define a circle xp-rangeError, yp, dpq+rangeError，左外
    t = 0:0.01:2*pi;
    a = (dpq-rangeError)*sin(t)+(xp+rangeError);
    b = (dpq-rangeError)*cos(t)+yp;
    in_circle22_1 = inpolygon(Xoq1,Yoq1,a,b);
    index22_1 = find(in_circle22_1==0);

    %Define a circle xp-rangeError, yp, dpq+rangeError，右外
    t = 0:0.01:2*pi;
    a = (dpq-rangeError)*sin(t)+(xp-rangeError);
    b = (dpq-rangeError)*cos(t)+yp;
    in_circle22_2 = inpolygon(Xoq1,Yoq1,a,b);
    index22_2 = find(in_circle22_2==0);

    %Define a circle xp+rangeError, yp, dpq-rangeError，左内
    t = 0:0.01:2*pi;
    a = (dpq+rangeError)*sin(t)+(xp-rangeError);
    b = (dpq+rangeError)*cos(t)+yp;
    in_circle21_1 = inpolygon(Xoq1,Yoq1,a,b);
    index21_1 = find(in_circle21_1>0);

    %Define a circle xp+rangeError, yp, dpq-rangeError，右内
    t = 0:0.01:2*pi;
    a = (dpq+rangeError)*sin(t)+(xp+rangeError);
    b = (dpq+rangeError)*cos(t)+yp;
    in_circle21_2 = inpolygon(Xoq1,Yoq1,a,b);
    index21_2 = find(in_circle21_2>0);
    
    %边缘交点
    Xoq1 = Xoq1(intersect(union(index22_1,index22_2),union(index21_1,index21_2)));
    Yoq1 = Yoq1(intersect(union(index22_1,index22_2),union(index21_1,index21_2)));
%     plot(Xoq1, Yoq1, 'b*');

    Xq1 = [Xpq1, Xoq1];    
    Yq1 = [Ypq1, Yoq1];    
%     plot(Xq, Xq, 'b*');

    %get the possible region boder of si
    %大半圆环
    xoi1 = -(doi+rangeError):mesh_density:(doi+rangeError);%out
    yoi11 = sqrt((doi+rangeError).^2 - xoi1.^2);
    yoi12 = -sqrt((doi+rangeError).^2 - xoi1.^2);
    %小半圆环
    xoi2 = -(doi-rangeError):mesh_density:(doi-rangeError);%in
    yoi21 = sqrt((doi-rangeError).^2 - xoi2.^2);
    yoi22 = -sqrt((doi-rangeError).^2 - xoi2.^2);
    %线段1
    xoi31 = -(doi+rangeError):mesh_density:-(doi-rangeError);
    yoi31 = zeros(1,length(xoi31));
    %线段2
    xoi32 = (doi-rangeError):mesh_density:(doi+rangeError);
    yoi32 = zeros(1,length(xoi32));
    %边缘线
    Xoi1 = [xoi1, xoi2, xoi31, xoi32];
    Yoi1 = [yoi11, yoi21, yoi31, yoi32];
    Xoi2 = [xoi1, xoi2, xoi31, xoi32];
    Yoi2 = [yoi12, yoi22, yoi31, yoi32];
%     plot(Xoi1, Yoi1, 'b*');
%     plot(Xoi2, Yoi2, 'b*');

    %左大半圆环
    xpi11 = (dop-rangeError)-(dpi+rangeError):mesh_density:dop;
    ypi111 = sqrt((dpi+rangeError).^2 - (xpi11-(dop-rangeError)).^2);
    ypi112 = -sqrt((dpi+rangeError).^2 - (xpi11-(dop-rangeError)).^2);
    %右大半圆环
    xpi12 = dop:mesh_density:(dop+rangeError)+(dpi+rangeError);
    ypi121 = sqrt((dpi+rangeError).^2 - (xpi12-((dop+rangeError))).^2);
    ypi122 = -sqrt((dpi+rangeError).^2 - (xpi12-((dop+rangeError))).^2);
    %左小半圆环
    xpi21 = (dop+rangeError)-(dpi-rangeError):mesh_density:dop;
    ypi211 = sqrt((dpi-rangeError).^2 - (xpi21-(dop+rangeError)).^2);
    ypi212 = -sqrt((dpi-rangeError).^2 - (xpi21-(dop+rangeError)).^2);
    %右小半圆环
    xpi22 = dop:mesh_density:(dop-rangeError)+(dpi-rangeError);
    ypi221 = sqrt((dpi-rangeError).^2 - (xpi22-(dop-rangeError)).^2);
    ypi222 = -sqrt((dpi-rangeError).^2 - (xpi22-(dop-rangeError)).^2);
    %线段1
    xpi31 = (dop-rangeError)-(dpi+rangeError):mesh_density:(dop+rangeError)-(dpi-rangeError);
    ypi31 = zeros(1,length(xpi31));
    %线段2
    xpi32 = (dop-rangeError)+(dpi-rangeError):mesh_density:(dop+rangeError)+(dpi+rangeError);
    ypi32 = zeros(1,length(xpi32));
    %线段2
    Xpi1 = [xpi11, xpi12, xpi21, xpi22, xpi31, xpi32];
    Ypi1 = [ypi111, ypi121, ypi211, ypi221, ypi31, ypi32];
    Xpi2 = [xpi11, xpi12, xpi21, xpi22, xpi31, xpi32];
    Ypi2 = [ypi112, ypi122, ypi212, ypi222, ypi31, ypi32];
%     plot(Xpi1, Ypi1, 'b*');
%     plot(Xpi2, Ypi2, 'b*');
    
    %Define a circle xo, yo, doi-rangeError，外
    t = 0:0.01:2*pi;
    a = (doi-rangeError)*sin(t)+xo;
    b = (doi-rangeError)*cos(t)+yo;
    in_circle12_1 = inpolygon(Xpi1,Ypi1,a,b);
    in_circle12_2 = inpolygon(Xpi2,Ypi2,a,b);
    index12_1 = find(in_circle12_1==0);
    index12_2 = find(in_circle12_2==0);
    
    %Define a circle xo, yo, doi+rangeError，内
    t = 0:0.01:2*pi;
    a = (doi+rangeError)*sin(t)+xo;
    b = (doi+rangeError)*cos(t)+yo;
    in_circle11_1 = inpolygon(Xpi1,Ypi1,a,b);
    in_circle11_2 = inpolygon(Xpi2,Ypi2,a,b);
    index11_1 = find(in_circle11_1>0);
    index11_2 = find(in_circle11_2>0);

    %边缘交点
    Xpi1 = Xpi1(intersect(index12_1, index11_1));%a part of region border
    Ypi1 = Ypi1(intersect(index12_1, index11_1));
    Xpi2 = Xpi2(intersect(index12_2, index11_2));%a part of region border
    Ypi2 = Ypi2(intersect(index12_2, index11_2));
%     plot(Xpi1, Ypi1, 'b*');
%     plot(Xpi2, Ypi2, 'b*');

    %Define a circle xp-rangeError, yp, dpi+rangeError，左外
    t = 0:0.01:2*pi;
    a = (dpi-rangeError)*sin(t)+(xp+rangeError);
    b = (dpi-rangeError)*cos(t)+yp;
    in_circle22_11 = inpolygon(Xoi1,Yoi1,a,b);
    in_circle22_12 = inpolygon(Xoi2,Yoi2,a,b);
    index22_11 = find(in_circle22_11==0);
    index22_12 = find(in_circle22_12==0);

    %Define a circle xp-rangeError, yp, dpi+rangeError，右外
    t = 0:0.01:2*pi;
    a = (dpi-rangeError)*sin(t)+(xp-rangeError);
    b = (dpi-rangeError)*cos(t)+yp;
    in_circle22_21 = inpolygon(Xoi1,Yoi1,a,b);
    in_circle22_22 = inpolygon(Xoi2,Yoi2,a,b);
    index22_21 = find(in_circle22_21==0);
    index22_22 = find(in_circle22_22==0);

    %Define a circle xp+rangeError, yp, dpi-rangeError，左内
    t = 0:0.01:2*pi;
    a = (dpi+rangeError)*sin(t)+(xp-rangeError);
    b = (dpi+rangeError)*cos(t)+yp;
    in_circle21_11 = inpolygon(Xoi1,Yoi1,a,b);
    in_circle21_12 = inpolygon(Xoi2,Yoi2,a,b);
    index21_11 = find(in_circle21_11>0);
    index21_12 = find(in_circle21_12>0);

    %Define a circle xp+rangeError, yp, dpi-rangeError，右内
    t = 0:0.01:2*pi;
    a = (dpi+rangeError)*sin(t)+(xp+rangeError);
    b = (dpi+rangeError)*cos(t)+yp;
    in_circle21_21 = inpolygon(Xoi1,Yoi1,a,b);
    in_circle21_22 = inpolygon(Xoi2,Yoi2,a,b);
    index21_21 = find(in_circle21_21>0);
    index21_22 = find(in_circle21_22>0);
    
    %边缘交点
    Xoi1 = Xoi1(intersect(union(index22_11,index22_21),union(index21_11,index21_21)));
    Yoi1 = Yoi1(intersect(union(index22_11,index22_21),union(index21_11,index21_21)));
    Xoi2 = Xoi2(intersect(union(index22_12,index22_22),union(index21_12,index21_22)));
    Yoi2 = Yoi2(intersect(union(index22_12,index22_22),union(index21_12,index21_22)));
%     plot(Xoi1, Yoi1, 'b*');
%     plot(Xoi2, Yoi2, 'b*');

    Xi1 = [Xpi1, Xoi1];    
    Yi1 = [Ypi1, Yoi1];    
    Xi2 = [Xpi2, Xoi2];    
    Yi2 = [Ypi2, Yoi2];    
%     plot(Xi1, Yi1, 'b*');
%     plot(Xi2, Yi2, 'b*');

    %判断两个区域之间距离的关系
    %Dqi1 = (Xq1-Xi1)^2 + (Yq1-Yi1)^2;
    %Dqi2 = (Xq1-Xi2)^2 + (Yq1-Yi2)^2;
    Dqi1 = [];
    for i = 1:length(Xq1),
        for j = 1:length(Xi1),
            dqi1 = sqrt((Xq1(i)-Xi1(j))^2 + (Yq1(i)-Yi1(j))^2);
            Dqi1 = [Dqi1; [dqi1]];
        end
    end
    Dqi2 = [];
    for i = 1:length(Xq1),
        for j = 1:length(Xi2),
            dqi2 = sqrt((Xq1(i)-Xi2(j))^2 + (Yq1(i)-Yi2(j))^2);
            Dqi2 = [Dqi2; [dqi2]];
        end
    end
    
    %if (Dqi1 < dqi-rangeError) && (Dqi1 > dqi+rangeError)... 
    %&& (Dqi2 > dqi-rangeError) && (Dqi2 < dqi+rangeError),
    if (length(find(Dqi1<dqi-rangeError)) == length(Dqi1) || length(find(Dqi1>dqi+rangeError)) == length(Dqi1))... %全不 All No
    && (~isempty(intersect(find(Dqi2>dqi-rangeError),find(Dqi2<dqi+rangeError)))),                 %部是 Part Yes
        interarea.area1 = 0;
        interarea.area2 = 1;
        return;
    %elseif (Dqi1 > dqi-rangeError) && (Dqi1 < dqi+rangeError)... 
    %    && (Dqi2 < dqi-rangeError) && (Dqi2 > dqi+rangeError),
    elseif (length(find(Dqi2<dqi-rangeError)) == length(Dqi2) || length(find(Dqi2>dqi+rangeError)) == length(Dqi2))... %全不 All No
        && (~isempty(intersect(find(Dqi1>dqi-rangeError),find(Dqi1<dqi+rangeError)))),                 %部是 Part Yes
        interarea.area1 = 1;
        interarea.area2 = 0;
        return;
    else
        interarea.area1 = -1;
        interarea.area2 = -1;
        return;
    end
end


% --- calculate the interarea by criterion2-bilateration
function interarea = GetCircleIntersecByCriterion2(xo, yo, xp, yp, xq, yq, dop, doq, dpq, doi, dpi, dr, rangeError, rangeVelocity)
disp(['xo:', num2str(xo),' yo:', num2str(yo),' xp:', num2str(xp),' yp:', num2str(yp),' xq:', num2str(xq),' yq:', num2str(yq),' dop:', num2str(dop),' doq:', num2str(doq),' dpq:', num2str(dpq),' doi:', num2str(doi),' dpi:', num2str(dpi),' rangeError:', num2str(rangeError)]);

%把速度和方向看成另一种误差因素参与运算
rangeError = rangeError + rangeVelocity;

if rangeError == 0,
    %the distance between so and sp
    r = sqrt((xp-xo)^2+(yp-yo)^2);%圆心距
    %如果两圆（oi与pi）相交
    if (doi+dpi>r) && abs(dpi-doi)<r,
        disp('相交1');
        %计算两圆的交点
        [xi,yi]=circcirc(xo, yo, doi, xp, yp, dpi);
        data = sortrows(double(vpa([xi;yi]', 5)),2)';
        xi1 = data(1,2); yi1 = data(2,2);
        xi2 = data(1,1); yi2 = data(2,1);
        %计算两个交点与sq的距离
        dqi1 = sqrt((xq-xi1)^2+(yq-yi1)^2);
        dqi2 = sqrt((xq-xi2)^2+(yq-yi2)^2);
        %确定哪一个交点是正确的
        if (dqi1 < dr) && (dqi2 > dr),
            interarea.area1 = 0;
            interarea.area2 = 1;
            return;
        %这个条件应该不会执行
        elseif (dqi1 > dr) && (dqi2 < dr),
            interarea.area1 = 1;
            interarea.area2 = 0;
            return;
        else
            interarea.area1 = -1;
            interarea.area2 = -1;
            return;
        end
    else
        disp('两圆不相交1');
        interarea.area1 = -1;
        interarea.area2 = -1;
        return;
    end
else
    mesh_density=rangeError/5;
    
    %get the possible region boder of sq
    %大半圆环
    xoq1 = -(doq+rangeError):mesh_density:(doq+rangeError);%out
    yoq1 = sqrt((doq+rangeError).^2 - xoq1.^2);
    %小半圆环
    xoq2 = -(doq-rangeError):mesh_density:(doq-rangeError);%in
    yoq2 = sqrt((doq-rangeError).^2 - xoq2.^2);
    %线段1
    xoq31 = -(doq+rangeError):mesh_density:-(doq-rangeError);
    yoq31 = zeros(1,length(xoq31));
    %线段2
    xoq32 = (doq-rangeError):mesh_density:(doq+rangeError);
    yoq32 = zeros(1,length(xoq32));
    %边缘线
    Xoq1 = [xoq1, xoq2, xoq31, xoq32];
    Yoq1 = [yoq1, yoq2, yoq31, yoq32];
%     plot(Xoq1, Yoq1, 'b*');

    %左大半圆环
    xpq11 = (dop-rangeError)-(dpq+rangeError):mesh_density:dop;
    ypq11 = sqrt((dpq+rangeError).^2 - (xpq11-(dop-rangeError)).^2);
    %右大半圆环
    xpq12 = dop:mesh_density:(dop+rangeError)+(dpq+rangeError);
    ypq12 = sqrt((dpq+rangeError).^2 - (xpq12-((dop+rangeError))).^2);
    %左小半圆环
    xpq21 = (dop+rangeError)-(dpq-rangeError):mesh_density:dop;
    ypq21 = sqrt((dpq-rangeError).^2 - (xpq21-(dop+rangeError)).^2);
    %右小半圆环
    xpq22 = dop:mesh_density:(dop-rangeError)+(dpq-rangeError);
    ypq22 = sqrt((dpq-rangeError).^2 - (xpq22-(dop-rangeError)).^2);
    %线段1
    xpq31 = (dop-rangeError)-(dpq+rangeError):mesh_density:(dop+rangeError)-(dpq-rangeError);
    ypq31 = zeros(1,length(xpq31));
    %线段2
    xpq32 = (dop-rangeError)+(dpq-rangeError):mesh_density:(dop+rangeError)+(dpq+rangeError);
    ypq32 = zeros(1,length(xpq32));
    %线段2
    Xpq1 = [xpq11, xpq12, xpq21, xpq22, xpq31, xpq32];
    Ypq1 = [ypq11, ypq12, ypq21, ypq22, ypq31, ypq32];
%     plot(Xpq1, Ypq1, 'b*');
    
    %Define a circle xo, yo, doq-rangeError，外
    t = 0:0.01:2*pi;
    a = (doq-rangeError)*sin(t)+xo;
    b = (doq-rangeError)*cos(t)+yo;
    in_circle12_1 = inpolygon(Xpq1,Ypq1,a,b);
    index12_1 = find(in_circle12_1==0);
    
    %Define a circle xo, yo, doq+rangeError，内
    t = 0:0.01:2*pi;
    a = (doq+rangeError)*sin(t)+xo;
    b = (doq+rangeError)*cos(t)+yo;
    in_circle11_1 = inpolygon(Xpq1,Ypq1,a,b);
    index11_1 = find(in_circle11_1>0);

    %边缘交点
    Xpq1 = Xpq1(intersect(index12_1, index11_1));%a part of region border
    Ypq1 = Ypq1(intersect(index12_1, index11_1));
%     plot(Xpq1, Ypq1, 'b*');

    %Define a circle xp-rangeError, yp, dpq+rangeError，左外
    t = 0:0.01:2*pi;
    a = (dpq-rangeError)*sin(t)+(xp+rangeError);
    b = (dpq-rangeError)*cos(t)+yp;
    in_circle22_1 = inpolygon(Xoq1,Yoq1,a,b);
    index22_1 = find(in_circle22_1==0);

    %Define a circle xp-rangeError, yp, dpq+rangeError，右外
    t = 0:0.01:2*pi;
    a = (dpq-rangeError)*sin(t)+(xp-rangeError);
    b = (dpq-rangeError)*cos(t)+yp;
    in_circle22_2 = inpolygon(Xoq1,Yoq1,a,b);
    index22_2 = find(in_circle22_2==0);

    %Define a circle xp+rangeError, yp, dpq-rangeError，左内
    t = 0:0.01:2*pi;
    a = (dpq+rangeError)*sin(t)+(xp-rangeError);
    b = (dpq+rangeError)*cos(t)+yp;
    in_circle21_1 = inpolygon(Xoq1,Yoq1,a,b);
    index21_1 = find(in_circle21_1>0);

    %Define a circle xp+rangeError, yp, dpq-rangeError，右内
    t = 0:0.01:2*pi;
    a = (dpq+rangeError)*sin(t)+(xp+rangeError);
    b = (dpq+rangeError)*cos(t)+yp;
    in_circle21_2 = inpolygon(Xoq1,Yoq1,a,b);
    index21_2 = find(in_circle21_2>0);
    
    %边缘交点
    Xoq1 = Xoq1(intersect(union(index22_1,index22_2),union(index21_1,index21_2)));
    Yoq1 = Yoq1(intersect(union(index22_1,index22_2),union(index21_1,index21_2)));
%     plot(Xoq1, Yoq1, 'b*');

    Xq1 = [Xpq1, Xoq1];    
    Yq1 = [Ypq1, Yoq1];    
%     plot(Xq, Xq, 'b*');

    %get the possible region boder of si
    %大半圆环
    xoi1 = -(doi+rangeError):mesh_density:(doi+rangeError);%out
    yoi11 = sqrt((doi+rangeError).^2 - xoi1.^2);
    yoi12 = -sqrt((doi+rangeError).^2 - xoi1.^2);
    %小半圆环
    xoi2 = -(doi-rangeError):mesh_density:(doi-rangeError);%in
    yoi21 = sqrt((doi-rangeError).^2 - xoi2.^2);
    yoi22 = -sqrt((doi-rangeError).^2 - xoi2.^2);
    %线段1
    xoi31 = -(doi+rangeError):mesh_density:-(doi-rangeError);
    yoi31 = zeros(1,length(xoi31));
    %线段2
    xoi32 = (doi-rangeError):mesh_density:(doi+rangeError);
    yoi32 = zeros(1,length(xoi32));
    %边缘线
    Xoi1 = [xoi1, xoi2, xoi31, xoi32];
    Yoi1 = [yoi11, yoi21, yoi31, yoi32];
    Xoi2 = [xoi1, xoi2, xoi31, xoi32];
    Yoi2 = [yoi12, yoi22, yoi31, yoi32];
%     plot(Xoi1, Yoi1, 'b*');
%     plot(Xoi2, Yoi2, 'b*');

    %左大半圆环
    xpi11 = (dop-rangeError)-(dpi+rangeError):mesh_density:dop;
    ypi111 = sqrt((dpi+rangeError).^2 - (xpi11-(dop-rangeError)).^2);
    ypi112 = -sqrt((dpi+rangeError).^2 - (xpi11-(dop-rangeError)).^2);
    %右大半圆环
    xpi12 = dop:mesh_density:(dop+rangeError)+(dpi+rangeError);
    ypi121 = sqrt((dpi+rangeError).^2 - (xpi12-((dop+rangeError))).^2);
    ypi122 = -sqrt((dpi+rangeError).^2 - (xpi12-((dop+rangeError))).^2);
    %左小半圆环
    xpi21 = (dop+rangeError)-(dpi-rangeError):mesh_density:dop;
    ypi211 = sqrt((dpi-rangeError).^2 - (xpi21-(dop+rangeError)).^2);
    ypi212 = -sqrt((dpi-rangeError).^2 - (xpi21-(dop+rangeError)).^2);
    %右小半圆环
    xpi22 = dop:mesh_density:(dop-rangeError)+(dpi-rangeError);
    ypi221 = sqrt((dpi-rangeError).^2 - (xpi22-(dop-rangeError)).^2);
    ypi222 = -sqrt((dpi-rangeError).^2 - (xpi22-(dop-rangeError)).^2);
    %线段1
    xpi31 = (dop-rangeError)-(dpi+rangeError):mesh_density:(dop+rangeError)-(dpi-rangeError);
    ypi31 = zeros(1,length(xpi31));
    %线段2
    xpi32 = (dop-rangeError)+(dpi-rangeError):mesh_density:(dop+rangeError)+(dpi+rangeError);
    ypi32 = zeros(1,length(xpi32));
    %线段2
    Xpi1 = [xpi11, xpi12, xpi21, xpi22, xpi31, xpi32];
    Ypi1 = [ypi111, ypi121, ypi211, ypi221, ypi31, ypi32];
    Xpi2 = [xpi11, xpi12, xpi21, xpi22, xpi31, xpi32];
    Ypi2 = [ypi112, ypi122, ypi212, ypi222, ypi31, ypi32];
%     plot(Xpi1, Ypi1, 'b*');
%     plot(Xpi2, Ypi2, 'b*');
    
    %Define a circle xo, yo, doi-rangeError，外
    t = 0:0.01:2*pi;
    a = (doi-rangeError)*sin(t)+xo;
    b = (doi-rangeError)*cos(t)+yo;
    in_circle12_1 = inpolygon(Xpi1,Ypi1,a,b);
    in_circle12_2 = inpolygon(Xpi2,Ypi2,a,b);
    index12_1 = find(in_circle12_1==0);
    index12_2 = find(in_circle12_2==0);
    
    %Define a circle xo, yo, doi+rangeError，内
    t = 0:0.01:2*pi;
    a = (doi+rangeError)*sin(t)+xo;
    b = (doi+rangeError)*cos(t)+yo;
    in_circle11_1 = inpolygon(Xpi1,Ypi1,a,b);
    in_circle11_2 = inpolygon(Xpi2,Ypi2,a,b);
    index11_1 = find(in_circle11_1>0);
    index11_2 = find(in_circle11_2>0);

    %边缘交点
    Xpi1 = Xpi1(intersect(index12_1, index11_1));%a part of region border
    Ypi1 = Ypi1(intersect(index12_1, index11_1));
    Xpi2 = Xpi2(intersect(index12_2, index11_2));%a part of region border
    Ypi2 = Ypi2(intersect(index12_2, index11_2));
%     plot(Xpi1, Ypi1, 'b*');
%     plot(Xpi2, Ypi2, 'b*');

    %Define a circle xp-rangeError, yp, dpi+rangeError，左外
    t = 0:0.01:2*pi;
    a = (dpi-rangeError)*sin(t)+(xp+rangeError);
    b = (dpi-rangeError)*cos(t)+yp;
    in_circle22_11 = inpolygon(Xoi1,Yoi1,a,b);
    in_circle22_12 = inpolygon(Xoi2,Yoi2,a,b);
    index22_11 = find(in_circle22_11==0);
    index22_12 = find(in_circle22_12==0);

    %Define a circle xp-rangeError, yp, dpi+rangeError，右外
    t = 0:0.01:2*pi;
    a = (dpi-rangeError)*sin(t)+(xp-rangeError);
    b = (dpi-rangeError)*cos(t)+yp;
    in_circle22_21 = inpolygon(Xoi1,Yoi1,a,b);
    in_circle22_22 = inpolygon(Xoi2,Yoi2,a,b);
    index22_21 = find(in_circle22_21==0);
    index22_22 = find(in_circle22_22==0);

    %Define a circle xp+rangeError, yp, dpi-rangeError，左内
    t = 0:0.01:2*pi;
    a = (dpi+rangeError)*sin(t)+(xp-rangeError);
    b = (dpi+rangeError)*cos(t)+yp;
    in_circle21_11 = inpolygon(Xoi1,Yoi1,a,b);
    in_circle21_12 = inpolygon(Xoi2,Yoi2,a,b);
    index21_11 = find(in_circle21_11>0);
    index21_12 = find(in_circle21_12>0);

    %Define a circle xp+rangeError, yp, dpi-rangeError，右内
    t = 0:0.01:2*pi;
    a = (dpi+rangeError)*sin(t)+(xp+rangeError);
    b = (dpi+rangeError)*cos(t)+yp;
    in_circle21_21 = inpolygon(Xoi1,Yoi1,a,b);
    in_circle21_22 = inpolygon(Xoi2,Yoi2,a,b);
    index21_21 = find(in_circle21_21>0);
    index21_22 = find(in_circle21_22>0);
    
    %边缘交点
    Xoi1 = Xoi1(intersect(union(index22_11,index22_21),union(index21_11,index21_21)));
    Yoi1 = Yoi1(intersect(union(index22_11,index22_21),union(index21_11,index21_21)));
    Xoi2 = Xoi2(intersect(union(index22_12,index22_22),union(index21_12,index21_22)));
    Yoi2 = Yoi2(intersect(union(index22_12,index22_22),union(index21_12,index21_22)));
%     plot(Xoi1, Yoi1, 'b*');
%     plot(Xoi2, Yoi2, 'b*');

    Xi1 = [Xpi1, Xoi1];    
    Yi1 = [Ypi1, Yoi1];    
    Xi2 = [Xpi2, Xoi2];    
    Yi2 = [Ypi2, Yoi2];    
%     plot(Xi1, Yi1, 'b*');
%     plot(Xi2, Yi2, 'b*');
    
    %判断两个区域之间距离的关系
    %Dqi1 = (Xq1-Xi1)^2 + (Yq1-Yi1)^2;
    %Dqi2 = (Xq1-Xi2)^2 + (Yq1-Yi2)^2;
    Dqi1 = [];
    for i = 1:length(Xq1),
        for j = 1:length(Xi1),
            dqi1 = sqrt((Xq1(i)-Xi1(j))^2 + (Yq1(i)-Yi1(j))^2);
            Dqi1 = [Dqi1; [dqi1]];
        end
    end
    Dqi2 = [];
    for i = 1:length(Xq1),
        for j = 1:length(Xi2),
            dqi2 = sqrt((Xq1(i)-Xi2(j))^2 + (Yq1(i)-Yi2(j))^2);
            Dqi2 = [Dqi2; [dqi2]];
        end
    end
    
    %if (Dqi1 < dr) && (Dqi2 > dr),
    if (length(find(Dqi1 < dr)) == length(Dqi1)) && (~isempty(find(Dqi2 > dr))),
        interarea.area1 = 0;
        interarea.area2 = 1;
        return;
    %这个条件应该不会执行
    %elseif (dqi2 < dr) && (dqi1 > dr),
    elseif (length(find(Dqi2 < dr)) == length(Dqi2)) && (~isempty(find(Dqi1 > dr))),
        interarea.area1 = 1;
        interarea.area2 = 0;
        return;
    else
        interarea.area1 = -1;
        interarea.area2 = -1;
        return;
    end
end


% --- bitri flip localization
function localCoor = biTriflipLocalize(vertex_o, vertex_p, vertex_q, vertex_i, dMatrix, mMatrix, distMatrix, tdistMatrix, measMatrix, rangeError, rangeVelocity)
% vertex_o, vertex_p, vertex_q, vertex_i, distMatrix, measMatrix.
%localCoor = struct('vertex_o',{},'vertex_o_x',{},'vertex_o_y',{},'vertex_p',{},'vertex_p_x',{},'vertex_p_y',{},'vertex_q',{},'vertex_q_x',{},'vertex_q_y',{},'nbVertexPosProb',{});
%calculate the position of vertex o,p,q

localCoor.nbVertexPosProb = -1;
localCoor.vertex_o = vertex_o;
localCoor.vertex_p = vertex_p;
localCoor.vertex_q = vertex_q;
localCoor.vertex_o_x = 0;
localCoor.vertex_o_y = 0;
localCoor.vertex_p_x = sqrt(dMatrix(vertex_o,vertex_p));
localCoor.vertex_p_y = 0;
localCoor.vertex_q_x = (dMatrix(vertex_o,vertex_q) + dMatrix(vertex_o,vertex_p) - dMatrix(vertex_p,vertex_q))/(2*sqrt(dMatrix(vertex_o,vertex_p)));
localCoor.vertex_q_y = sqrt(dMatrix(vertex_o,vertex_q) - (localCoor.vertex_q_x).^2);

%calculate the probabile position of vertex i
localCoor.vertex_i = vertex_i;
localCoor.vertex_i_x  = (dMatrix(vertex_o,vertex_i) + dMatrix(vertex_o,vertex_p) - dMatrix(vertex_p,vertex_i))/(2*sqrt(dMatrix(vertex_o,vertex_p)));
vertex_i_y1 = +sqrt(dMatrix(vertex_o,vertex_i) - (localCoor.vertex_i_x).^2);
vertex_i_y2 = -sqrt(dMatrix(vertex_o,vertex_i) - (localCoor.vertex_i_x).^2);

if mMatrix(vertex_i,vertex_q) == 1, %trilateration
   interarea = GetCircleIntersecByCriterion1(localCoor.vertex_o_x, localCoor.vertex_o_y, localCoor.vertex_p_x, localCoor.vertex_p_y, localCoor.vertex_q_x, localCoor.vertex_q_y,...
                                             sqrt(dMatrix(vertex_o,vertex_p)), sqrt(dMatrix(vertex_o,vertex_q)), sqrt(dMatrix(vertex_p,vertex_q)),...
                                             sqrt(dMatrix(vertex_o,vertex_i)), sqrt(dMatrix(vertex_p,vertex_i)), sqrt(dMatrix(vertex_q,vertex_i)),...
                                             rangeError, rangeVelocity);
   if (interarea.area1 > 0) && (interarea.area2 == 0),
      %display(['(i)',num2str(vertex_i),' (', num2str(dMatrix(vertex_i,vertex_q)),'=',num2str(dSqr1),')',' and (', num2str(dMatrix(vertex_i,vertex_q)),'~=',num2str(dSqr2),')']);
            
      localCoor.vertex_i_y = vertex_i_y1;
      localCoor.nbVertexPosProb = 1;
   elseif (interarea.area1 == 0) && (interarea.area2 > 0),
      %display(['(i)',num2str(vertex_i),' (', num2str(dMatrix(vertex_i,vertex_q)),'=',num2str(dSqr2),')',' and (', num2str(dMatrix(vertex_i,vertex_q)),'~=',num2str(dSqr1),')']);

      localCoor.vertex_i_y = vertex_i_y2;
      localCoor.nbVertexPosProb = 1;
   else
      localCoor.nbVertexPosProb = 0;
   end 
else %考虑通讯range的两点定位
   interarea = GetCircleIntersecByCriterion2(localCoor.vertex_o_x, localCoor.vertex_o_y, localCoor.vertex_p_x, localCoor.vertex_p_y, localCoor.vertex_q_x, localCoor.vertex_q_y,...
                                             sqrt(dMatrix(vertex_o,vertex_p)), sqrt(dMatrix(vertex_o,vertex_q)), sqrt(dMatrix(vertex_p,vertex_q)),...
                                             sqrt(dMatrix(vertex_o,vertex_i)), sqrt(dMatrix(vertex_p,vertex_i)),...
                                             min([max(sqrt(distMatrix(vertex_q,:))),max(sqrt(distMatrix(vertex_i,:)))]),...
                                             rangeError, rangeVelocity);
   if (interarea.area1 > 0) && (interarea.area2 == 0),
      %display(['(i)',num2str(vertex_i),' (', num2str(dSqr1),'<=',num2str(maxRangeRadiusSqr),')',' and (', num2str(dSqr2),'>',num2str(maxRangeRadiusSqr),')']);

      localCoor.vertex_i_y = vertex_i_y1;
      localCoor.nbVertexPosProb = 1;
   elseif (interarea.area1 == 0) && (interarea.area2 > 0),
      %display(['(i)',num2str(vertex_i),' (', num2str(dSqr2),'<=',num2str(maxRangeRadiusSqr),')',' and (', num2str(dSqr1),'>',num2str(maxRangeRadiusSqr),')']);

      localCoor.vertex_i_y = vertex_i_y2;
      localCoor.nbVertexPosProb = 1;
   else
      localCoor.nbVertexPosProb = 0;
   end
end


% --- biTrilateration localization
function [xl, yl, vertexSeqLocl] = biTrilateration(networkidx,graphidx,total,x0,y0,vertex_fo,vertex_fp,vertex_fq,rangeRadius,rangeError,rangeVelocity,rangeShield,distMatrix,tdistMatrix,measMatrix,degreeSeq,vertexSeq,nbDegreeSeq,nbVertexSeq);

t1 = clock;
display(['Arithmetic begin time:', datestr(clock,0)]);

plot(x0([vertex_fo,vertex_fp,vertex_fq]), y0([vertex_fo,vertex_fp,vertex_fq]), 'b*');

vertexAll = 1:total;
vertexSeqLocl = [vertex_fo,vertex_fp,vertex_fq];
vertexSequnLocl = setdiff(vertexAll, vertexSeqLocl);

vertex_o = vertex_fo; %o
vertex_p = vertex_fp; %p
vertex_q = vertex_fq; %q
         
xl(vertex_o) = x0(vertex_o);
yl(vertex_o) = y0(vertex_o);
xl(vertex_p) = x0(vertex_p);
yl(vertex_p) = y0(vertex_p);
xl(vertex_q) = x0(vertex_q);
yl(vertex_q) = y0(vertex_q);

newTriangle = cell(1,1);
newTriangle{1} = [newTriangle{1},[vertex_o;vertex_p;vertex_q]];

while(1)
    vertexSeqLocl1 = vertexSeqLocl;
    
    for i = 1:length(vertexSequnLocl),
        vertex_i = vertexSequnLocl(i);
        
        iVertexSeq = nbVertexSeq(vertex_i,:); %neighbor of nonlocalized
        tag = 0;
        for j = 1:length(newTriangle),
            vertexSame = intersect(newTriangle{j}',iVertexSeq);

            %the number of the same vertex >= 3
            if length(vertexSame) >= 3,
                display(['i:', num2str(vertex_i), ':vertexSame:', mat2str(vertexSame)]);

                flag=0;
                for oi = 1:size(vertexSame,2);
                    vertex_o = vertexSame(oi);
                    if vertex_o == vertex_i,
                        continue;
                    end
                
                    for pi = (oi+1):size(vertexSame,2);
                        vertex_p = vertexSame(pi);
                        if vertex_p == vertex_i || vertex_p == vertex_o,
                            continue;
                        end
                    
                        for qi = (pi+1):size(vertexSame,2);
                            vertex_q = vertexSame(qi);
                            if vertex_q == vertex_i || vertex_q == vertex_o || vertex_q == vertex_p,
                                continue;
                            end

                            if ~((measMatrix(vertex_o,vertex_p) == 1) && (measMatrix(vertex_o,vertex_q) == 1) && (measMatrix(vertex_p,vertex_q) == 1)...
                            && IsTriangle(distMatrix(vertex_o,vertex_p),distMatrix(vertex_o,vertex_q),distMatrix(vertex_p,vertex_q))),
                                display(['(o,p,q)',num2str(vertex_o), ',',num2str(vertex_p), ',',num2str(vertex_q), ' is not triangle']);
                                continue;
                            end
                            if ~((measMatrix(vertex_o,vertex_i) == 1) && (measMatrix(vertex_p,vertex_i) == 1) && (measMatrix(vertex_o,vertex_p) == 1)...
                            && IsTriangle(distMatrix(vertex_o,vertex_i),distMatrix(vertex_p,vertex_i),distMatrix(vertex_o,vertex_p))),
                                display(['(o,p,i)',num2str(vertex_o), ',',num2str(vertex_p), ',',num2str(vertex_i), ' is not triangle']);
                                continue;
                            end
                        
                            lc = biTriflipLocalize(vertex_o, vertex_p, vertex_q, vertex_i, distMatrix, measMatrix, distMatrix, tdistMatrix, measMatrix, rangeError*rangeShield, rangeVelocity);
                            if lc.nbVertexPosProb == 1,
                                display(['new triangle(o,p,q):', num2str(vertex_o), ',',num2str(vertex_p), ',',num2str(vertex_q)]);
                                display(['new vertex(i)',num2str(vertex_i)]);

                                vertexSeqLocl1 = union(vertexSeqLocl1, vertex_i);

                                M = [lc.vertex_o_x,lc.vertex_o_y,1; lc.vertex_p_x,lc.vertex_p_y,1; lc.vertex_q_x,lc.vertex_q_y,1]\...
                                    [xl(vertex_o),yl(vertex_o),1; xl(vertex_p),yl(vertex_p),1; xl(vertex_q),yl(vertex_q),1];
                                z = [lc.vertex_i_x,lc.vertex_i_y,1] * M;
                                xl(vertex_i) = z(1);
                                yl(vertex_i) = z(2);
                                    
                                newTriangle{length(newTriangle)+1} = [vertex_o;vertex_p;vertex_i];

                                %modify the distance of the localized subgraph
                                for li = 1:(total-1),
                                    for lj = (li+1):total,
                                        if ~isempty(find(vertexSeqLocl1 == li)) && ~isempty(find(vertexSeqLocl1 == lj)) && (measMatrix(li,lj) == 1),
                                            distMatrix(li,lj) = distVertexes(xl(li),yl(li),xl(lj),yl(lj));
                                            distMatrix(lj,li) = distMatrix(li,lj);
                                        end
                                    end
                                end
                                
                                plot(xl(vertex_i), yl(vertex_i), 'g*');
                                xp = x0(vertex_i):(xl(vertex_i)-x0(vertex_i))/10:xl(vertex_i);
                                yp = ((y0(vertex_i) - yl(vertex_i))/(x0(vertex_i) - xl(vertex_i)))*(xp - xl(vertex_i)) + yl(vertex_i);
                                plot(xp, yp, 'k');

                                flag =1;
                                tag = 1;
                                break;
                            end
                        end
                        if flag == 1,
                            break;
                        end;
                    end
                    if flag == 1,
                        break;
                    end;
                end
            
                if flag == 0,
                    display(['i:', num2str(vertex_i), ':vertexSame:', mat2str(vertexSame)]);

                    flag=0;
                    for oi = 1:size(vertexSame,2);
                        vertex_o = vertexSame(oi);
                        if vertex_o == vertex_i,
                            continue;
                        end
                
                        for pi = (oi+1):size(vertexSame,2);
                            vertex_p = vertexSame(pi);
                            if vertex_p == vertex_i || vertex_p == vertex_o,
                                continue;
                            end
                    
                            for qi = 1:size(vertexSeqLocl1,2);
                                vertex_q = vertexSeqLocl1(qi);
                                if (vertex_q == vertex_o)||(vertex_q == vertex_p)||(vertex_q == vertex_i),
                                    continue;
                                end
                       
                                if ~((measMatrix(vertex_o,vertex_p) == 1) && (measMatrix(vertex_o,vertex_q) == 1) && (measMatrix(vertex_p,vertex_q) == 1)...
                                && IsTriangle(distMatrix(vertex_o,vertex_p),distMatrix(vertex_o,vertex_q),distMatrix(vertex_p,vertex_q))),
                                    display(['(o,p,q)',num2str(vertex_o), ',',num2str(vertex_p), ',',num2str(vertex_q), ' is not triangle']);
                                    continue;
                                end
                                if ~((measMatrix(vertex_o,vertex_i) == 1) && (measMatrix(vertex_p,vertex_i) == 1) && (measMatrix(vertex_o,vertex_p) == 1)...
                                && IsTriangle(distMatrix(vertex_o,vertex_i),distMatrix(vertex_p,vertex_i),distMatrix(vertex_o,vertex_p))),
                                    display(['(o,p,i)',num2str(vertex_o), ',',num2str(vertex_p), ',',num2str(vertex_i), ' is not triangle']);
                                    continue;
                                end
                        
                                lc = biTriflipLocalize(vertex_o, vertex_p, vertex_q, vertex_i, distMatrix, measMatrix, distMatrix, tdistMatrix, measMatrix, rangeError*rangeShield, rangeVelocity);
                                if lc.nbVertexPosProb == 1,
                                    display(['new triangle(o,p,q):', num2str(vertex_o), ',',num2str(vertex_p), ',',num2str(vertex_q)]);
                                    display(['new vertex(i)',num2str(vertex_i)]);
                                    
                                    vertexSeqLocl1 = union(vertexSeqLocl1, vertex_i);

                                    M = [lc.vertex_o_x,lc.vertex_o_y,1; lc.vertex_p_x,lc.vertex_p_y,1; lc.vertex_q_x,lc.vertex_q_y,1]\...
                                        [xl(vertex_o),yl(vertex_o),1; xl(vertex_p),yl(vertex_p),1; xl(vertex_q),yl(vertex_q),1];
                                    z = [lc.vertex_i_x,lc.vertex_i_y,1] * M;
                                    xl(vertex_i) = z(1);
                                    yl(vertex_i) = z(2);

                                    newTriangle{length(newTriangle)+1} = [vertex_o;vertex_p;vertex_i];
                                    
                                    %modify the distance of the localized subgraph
                                    for li = 1:(total-1),
                                        for lj = (li+1):total,
                                            if ~isempty(find(vertexSeqLocl1 == li)) && ~isempty(find(vertexSeqLocl1 == lj)) && (measMatrix(li,lj) == 1),
                                                distMatrix(li,lj) = distVertexes(xl(li),yl(li),xl(lj),yl(lj));
                                                distMatrix(lj,li) = distMatrix(li,lj);
                                            end
                                        end
                                    end
                                
                                    plot(xl(vertex_i), yl(vertex_i), 'g*');
                                    xp = x0(vertex_i):(xl(vertex_i)-x0(vertex_i))/10:xl(vertex_i);
                                    yp = ((y0(vertex_i) - yl(vertex_i))/(x0(vertex_i) - xl(vertex_i)))*(xp - xl(vertex_i)) + yl(vertex_i);
                                    plot(xp, yp, 'k');

                                    flag =1;
                                    tag = 1;
                                    break;
                                end
                            end
                            if flag == 1,
                                break;
                            end;
                        end
                        if flag == 1,
                            break;
                        end;
                    end
                end
            elseif length(vertexSame) == 2,
                display(['i:', num2str(vertex_i), ':vertexSame:', mat2str(vertexSame)]);

                flag=0;
                for oi = 1:size(vertexSame,2);
                    vertex_o = vertexSame(oi);
                    if vertex_o == vertex_i,
                        continue;
                    end
                
                    for pi = (oi+1):size(vertexSame,2);
                        vertex_p = vertexSame(pi);
                        if vertex_p == vertex_i || vertex_p == vertex_o,
                            continue;
                        end
                    
                        for qi = 1:size(vertexSeqLocl1,2);
                            vertex_q = vertexSeqLocl1(qi);
                            if vertex_q == vertex_i || vertex_q == vertex_o || vertex_q == vertex_p,
                                continue;
                            end
                        
                            if ~((measMatrix(vertex_o,vertex_p) == 1) && (measMatrix(vertex_o,vertex_q) == 1) && (measMatrix(vertex_p,vertex_q) == 1)...
                            && IsTriangle(distMatrix(vertex_o,vertex_p),distMatrix(vertex_o,vertex_q),distMatrix(vertex_p,vertex_q))),
                                display(['(o,p,q)',num2str(vertex_o), ',',num2str(vertex_p), ',',num2str(vertex_q), ' is not triangle']);
                                continue;
                            end
                            if ~((measMatrix(vertex_o,vertex_i) == 1) && (measMatrix(vertex_p,vertex_i) == 1) && (measMatrix(vertex_o,vertex_p) == 1)...
                            && IsTriangle(distMatrix(vertex_o,vertex_i),distMatrix(vertex_p,vertex_i),distMatrix(vertex_o,vertex_p))),
                                display(['(o,p,i)',num2str(vertex_o), ',',num2str(vertex_p), ',',num2str(vertex_i), ' is not triangle']);
                                continue;
                            end
                        
                            lc = biTriflipLocalize(vertex_o, vertex_p, vertex_q, vertex_i, distMatrix, measMatrix, distMatrix, tdistMatrix, measMatrix, rangeError*rangeShield, rangeVelocity);
                            if lc.nbVertexPosProb == 1,
                                display(['new triangle(o,p,q):', num2str(vertex_o), ',',num2str(vertex_p), ',',num2str(vertex_q)]);
                                display(['new vertex(i)',num2str(vertex_i)]);

                                vertexSeqLocl1 = union(vertexSeqLocl1, vertex_i);
                                
                                M = [lc.vertex_o_x,lc.vertex_o_y,1; lc.vertex_p_x,lc.vertex_p_y,1; lc.vertex_q_x,lc.vertex_q_y,1]\[xl(vertex_o),yl(vertex_o),1; xl(vertex_p),yl(vertex_p),1; xl(vertex_q),yl(vertex_q),1];
                                z = [lc.vertex_i_x,lc.vertex_i_y,1] * M;
                                xl(vertex_i) = z(1);
                                yl(vertex_i) = z(2);

                                newTriangle{length(newTriangle)+1} = [vertex_o;vertex_p;vertex_i];
                                
                                %modify the distance of the localized subgraph
                                for li = 1:(total-1),
                                    for lj = (li+1):total,
                                        if ~isempty(find(vertexSeqLocl1 == li)) && ~isempty(find(vertexSeqLocl1 == lj)) && (measMatrix(li,lj) == 1),
                                            distMatrix(li,lj) = distVertexes(xl(li),yl(li),xl(lj),yl(lj));
                                            distMatrix(lj,li) = distMatrix(li,lj);
                                        end
                                    end
                                end
                                
                                plot(xl(vertex_i), yl(vertex_i), 'g*');
                                xp = x0(vertex_i):(xl(vertex_i)-x0(vertex_i))/10:xl(vertex_i);
                                yp = ((y0(vertex_i) - yl(vertex_i))/(x0(vertex_i) - xl(vertex_i)))*(xp - xl(vertex_i)) + yl(vertex_i);
                                plot(xp, yp, 'k');

                                flag =1;
                                tag = 1;
                                break;
                            end
                        end
                        if flag == 1,
                            break;
                        end;
                    end
                    if flag == 1,
                        break;
                    end;
                end
            end
        
            if tag == 1,
                break;
            end;
        end
    end
    
    if isequal(vertexSeqLocl,vertexSeqLocl1),
        break;
    end
    
    vertexSeqLocl = vertexSeqLocl1;
    vertexSequnLocl = setdiff(vertexAll, vertexSeqLocl);
end

error = 0;
for i = 1:total,
   if ~isempty(find(vertexSeqLocl == i));
      error = error + sqrt((x0(i) - xl(i)).^2 + (y0(i) - yl(i)).^2);
   end
end

plot(x0(vertexSeqLocl), y0(vertexSeqLocl), 'r*');
plot(xl(vertexSeqLocl), yl(vertexSeqLocl), 'g*');

for i = 1:total,
   if length(find(vertexSeqLocl == i)) == 1;
     xp = x0(i):(xl(i)-x0(i))/10:xl(i);
     yp = ((y0(i) - yl(i))/(x0(i) - xl(i)))*(xp - xl(i)) + yl(i);
     plot(xp, yp, 'k');
   end
end
plot(x0([vertex_fo,vertex_fp,vertex_fq]), y0([vertex_fo,vertex_fp,vertex_fq]), 'b*');

t2 = clock;

disp('*********************************************************************');
display(['vertex total:', num2str(total)]);
display(['range radius:', num2str(rangeRadius)]);
display(['range error:', num2str(rangeError)]);
display(['range velocity:', num2str(rangeVelocity)]);
display(['range shield:', mat2str(rangeShield)]);
display(['degree seq:', mat2str(degreeSeq)]);
display(['vertex seq:', mat2str(vertexSeq)]);
display(['average vertex degree:', num2str(mean(degreeSeq(:)))]);
display(['located vertex seq:', mat2str(vertexSeqLocl)]);
display(['working time:', num2str(etime(t2,t1))]);
display(['x0:', mat2str(x0)]);
display(['y0:', mat2str(y0)]);
% display(['x:', mat2str(x)]);
% display(['y:', mat2str(y)]);
display(['xl:', mat2str(xl(vertexSeqLocl))]);
display(['yl:', mat2str(yl(vertexSeqLocl))]);
display(['located vertex count:', num2str(size(vertexSeqLocl,2))]);
display(['located vertex rate:', num2str(size(vertexSeqLocl,2)*100/total), '%']);
display(['located error:', num2str(error/size(vertexSeqLocl,2),5)]);
disp('*********************************************************************');

display(['Arithmetic over time:', datestr(clock,0)]);

fid=fopen('output.txt','a');
fprintf(fid,'%6.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.5f  \n',networkidx,3.1,graphidx,total,rangeRadius,mean(degreeSeq(:)),size(vertexSeqLocl,2),size(vertexSeqLocl,2)*100/total,rangeError,rangeVelocity,rangeShield,error/size(vertexSeqLocl,2));
fclose(fid);


% --- Cut Vertex Sub procedure.
function nc = calcCutVertexSet(G)
% G adjacency matrix.
% nc set of cut vertices.
% Example 4.8
% G=[0,1,0,0;
%    1,0,0,1;
%    0,0,0,1;
%    0,1,1,0];
%
nc = [];
time_step = 0;
sta = [];              % Stack is initially empty.
N = size(G, 1);
p = zeros(1, N);
f = zeros(1, N);
l = zeros(1, N);
g = zeros(1, N);
h = zeros(1, N);
f_s = zeros(1, N);

time_step = time_step + 1;
f(1) = time_step;
sta(1) = 1;                        % Current vertex is top of the stack.
top = length(sta);

while top ~= 0                     % Execute while stack is not empty.
    time_step = time_step + 1;
    a = intersect(find(G(sta(top),:) == 1), find(f == 0)); 
    % Find the unlabelled neighbous of the current vertex.
    if ~isempty(a)
        f(a(1)) = time_step;       % Add the neighbor to the top of stack.
        p(a(1)) = sta(top);
        top =top + 1;
        sta(top) = a(1);
    else                           % Remove the current vertex from stack.
        l(sta(top)) = time_step;
        T = union(find(p == sta(top)), p(sta(top)));
        B = setdiff(find(G(sta(top),:) == 1), T);
        T = setdiff(T, p(sta(top)));
        if isempty(T) && isempty(B)
            f_s(sta(top)) = f(sta(top));
        elseif isempty(T) && ~isempty(B)
            g(sta(top)) = min(f(B));
            f_s(sta(top)) = min([f(sta(top)), g(sta(top))]);
        elseif ~isempty(T) && isempty(B)
            h(sta(top)) = min(f_s(T));
            f_s(sta(top)) = min([f(sta(top)), h(sta(top))]);
        else
            g(sta(top)) = min(f(B));
            h(sta(top)) = min(f_s(T));
            f_s(sta(top)) = min([f(sta(top)), g(sta(top)), h(sta(top))]);
        end
        sta(top) = [];
        top = top - 1;
    end
end

if length(find(p == 1)) >= 2
    nc = 1;
end
for m = 2:N
    a = find(p == m);
    if ~isempty(a)
        for n = 1:length(a)
            if f(m) <= f_s(a(n))
                nc = union(nc, m);
            end
        end
    end
end


function edtVertexTotal_Callback(hObject, eventdata, handles)
% hObject    handle to edtVertexTotal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtVertexTotal as text
%        str2double(get(hObject,'String')) returns contents of edtVertexTotal as a double


% --- Executes during object creation, after setting all properties.
function edtVertexTotal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtVertexTotal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edtRangeRadius_Callback(hObject, eventdata, handles)
% hObject    handle to edtRangeRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtRangeRadius as text
%        str2double(get(hObject,'String')) returns contents of edtRangeRadius as a double


% --- Executes during object creation, after setting all properties.
function edtRangeRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtRangeRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edtRangeError_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtRangeError (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edtRangeError_Callback(hObject, eventdata, handles)
% hObject    handle to edtRangeError (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtRangeError as text
%        str2double(get(hObject,'String')) returns contents of edtRangeError as a double


% --- Executes during object creation, after setting all properties.
function edtRangeShield_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtRangeShield (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edtRangeShield_Callback(hObject, eventdata, handles)
% hObject    handle to edtRangeShield (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtRangeShield as text
%        str2double(get(hObject,'String')) returns contents of edtRangeShield as a double


% --- Executes on button press in btnCreateVertex.
function btnCreateVertex_Callback(hObject, eventdata, handles)
% hObject    handle to btnCreateVertex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global total;
global x0;
global y0;

clc; %clear screen
cla; %clear axes

total = str2double(get(handles.edtVertexTotal,'String'));

if get(handles.cbReadmat,'Value') ~= 1,
   rng('shuffle', 'v5uniform');
   %x0 = [unifrnd(0,total,total,1)'];
   x0 = [unifrnd(0,100,total,1)'];
   rng('shuffle', 'v5uniform');
   %y0 = [unifrnd(0,total,total,1)'];
   y0 = [unifrnd(0,100,total,1)'];
end

display(['total:', num2str(total)]);
display(['x0:', mat2str(x0)]);
display(['y0:', mat2str(y0)]);


% --- Executes on button press in btnCreateError.
function btnCreateError_Callback(hObject, eventdata, handles)
% hObject    handle to btnCreateError (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global total;
global x0;
global y0;
global xe;
global ye;
global x;
global y;
global rangeRadius;
global rangeError;
global rangeShield;

rangeRadius = str2double(get(handles.edtRangeRadius,'String'));
rangeError = str2double(get(handles.edtRangeError,'String'));
rangeShield = str2double(get(handles.edtRangeShield,'String'));

rng('shuffle', 'v5uniform');
errorAngle = unifrnd(0,2*pi,total,1);
rng('shuffle', 'v5normal');
errorDist = sqrt(0.5*rangeError)*randn(1,total);
errorPosition = [errorDist' errorDist'].*[cos(errorAngle) sin(errorAngle)];
xe = [errorPosition(:,1)'];
ye = [errorPosition(:,2)'];

x = x0 + xe;
y = y0 + ye;

display(['rangeError:', mat2str(rangeError)]);
display(['rangeShield:', mat2str(rangeShield)]);
display(['x0:', mat2str(x0)]);
display(['y0:', mat2str(y0)]);
display(['xe:', mat2str(xe)]);
display(['ye:', mat2str(ye)]);
display(['x:', mat2str(x)]);
display(['y:', mat2str(y)]);

set(handles.tbSeq, 'data', [x; y]');

axis([0 100 0 100]);
axis manual;
grid off; % grid on;
hold on; % hold off

title('the graph of random node（x，y）');
xlabel('X-axis');
ylabel('Y-axis');
plot(x, y, 'k*');


% --- Executes on button press in cbSavemat.
function cbSavemat_Callback(hObject, eventdata, handles)
% hObject    handle to cbSavemat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbSavemat

global total;
global x0;
global y0;
global rangeRadius;
global rangeError;
global rangeShield;

rangeRadius = str2double(get(handles.edtRangeRadius,'String'));
rangeError = str2double(get(handles.edtRangeError,'String'));
rangeShield = str2double(get(handles.edtRangeShield,'String'));

if get(handles.cbSavemat,'Value') == 1,
   save .\mat\total total;
   save .\mat\x0 x0;
   save .\mat\y0 y0;
   save .\mat\rangeRadius rangeRadius;
   save .\mat\rangeError rangeError;
   save .\mat\rangeShield rangeShield;
end


% --- Executes on button press in cbReadmat.
function cbReadmat_Callback(hObject, eventdata, handles)
% hObject    handle to cbReadmat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbReadmat

global total;
global x0;
global y0;
global rangeRadius;
global rangeError;
global rangeShield;

if get(handles.cbReadmat,'Value') == 1,
   load .\mat\total total;
   load .\mat\x0 x0;
   load .\mat\y0 y0;
   load .\mat\rangeRadius rangeRadius;
   load .\mat\rangeError rangeError;
   load .\mat\rangeShield rangeShield;

   set(handles.edtVertexTotal,'string', total);
   set(handles.tbSeq, 'data', [x0; y0]');
   set(handles.edtRangeRadius,'string', rangeRadius);
   set(handles.edtRangeError,'string', rangeError);
   set(handles.edtRangeShield,'string', rangeShield);
end


% --- Executes on button press in rbBlack.
function rbBlack_Callback(hObject, eventdata, handles)
% hObject    handle to rbBlack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbBlack


% --- Executes on button press in rbRed.
function rbRed_Callback(hObject, eventdata, handles)
% hObject    handle to rbRed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbRed


% --- Executes on button press in rbBlue.
function rbBlue_Callback(hObject, eventdata, handles)
% hObject    handle to rbBlue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbBlue


% --- Executes on button press in cbShowSingleVertex.
function cbShowSingleVertex_Callback(hObject, eventdata, handles)
% hObject    handle to cbShowSingleVertex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbShowSingleVertex


% --- Executes on button press in cbShowSingleLink.
function cbShowSingleLink_Callback(hObject, eventdata, handles)
% hObject    handle to cbShowSingleLink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbShowSingleLink


% --- Executes on button press in cbShowSingleRadius.
function cbShowSingleRadius_Callback(hObject, eventdata, handles)
% hObject    handle to cbShowSingleRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbShowSingleRadius


% --- Executes on button press in cbShowSingleCoor.
function cbShowSingleCoor_Callback(hObject, eventdata, handles)
% hObject    handle to cbShowSingleCoor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbShowSingleCoor


% --- Executes on button press in cbShowRotation.
function cbShowRotation_Callback(hObject, eventdata, handles)
% hObject    handle to cbShowRotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbShowRotation


function edtVertexSeq_Callback(hObject, eventdata, handles)
% hObject    handle to edtVertexSeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtVertexSeq as text
%        str2double(get(hObject,'String')) returns contents of edtVertexSeq as a double


% --- Executes during object creation, after setting all properties.
function edtVertexSeq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtVertexSeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.


% --- Executes on button press in btnShowVertex.
function btnShowVertex_Callback(hObject, eventdata, handles)
% hObject    handle to btnShowVertex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global x0;
global y0;

axis([0 100 0 100]);
axis manual;
grid off; % grid on;
hold on; % hold off

title('the graph of random node（x，y）');
xlabel('X-axis');
ylabel('Y-axis');

if get(handles.cbShowSingleVertex,'Value') == 1,
    k = str2double(get(handles.edtVertexSeq,'String'));
 
    if get(handles.rbRed,'Value') == 1,
      plot(x0(k), y0(k), 'r*');
    elseif get(handles.rbBlue,'Value') == 1,
      plot(x0(k), y0(k), 'b*');
    elseif get(handles.rbWhite,'Value') == 1,
      plot(x0(k), y0(k), 'w*');
    else
      plot(x0(k), y0(k), 'k*');
    end
else
    if get(handles.rbRed,'Value') == 1,
      plot(x0, y0, 'r*');
    elseif get(handles.rbBlue,'Value') == 1,
      plot(x0, y0, 'b*');
    elseif get(handles.rbWhite,'Value') == 1,
      plot(x0, y0, 'w*');
    else
      plot(x0, y0, 'k*');
    end
end


% --- Executes on button press in btnShowLink.
function btnShowLink_Callback(hObject, eventdata, handles)
% hObject    handle to btnShowLink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global total;
global x0;
global y0;
global x;
global y;
global rangeRadius;

rangeRadius = str2double(get(handles.edtRangeRadius,'String'));
rangeRadiusSqr = rangeRadius^2;

if get(handles.cbShowSingleLink,'Value') == 1,
    k = str2double(get(handles.edtVertexSeq,'String'));
 
    for el = 1:total,
      if el == k,
         continue;
      end
      
      dSqr = distVertexes(x0(k),y0(k),x0(el),y0(el));
      %dSqr = distVertexes(x(k),y(k),x(el),y(el));
      if (dSqr <= rangeRadiusSqr) == 1,
          % 根据节点通讯半径，划线
          %plot(x(k), y(k), x(el), y(el));
          xl = x(k):(x(el)-x(k))/10:x(el);
          yl = ((y(k) - y(el))/(x(k) - x(el)))*(xl - x(el)) + y(el);
          if get(handles.rbRed,'Value') == 1,
             plot(xl, yl, 'r');
          elseif get(handles.rbBlue,'Value') == 1,
             plot(xl, yl, 'b');
          elseif get(handles.rbWhite,'Value') == 1,
             plot(xl, yl, 'w');
          else
             plot(xl, yl, 'k');
          end    
      end
   end
else
   for k = 1:(total-1),
      for el = (k+1):total,
         dSqr = distVertexes(x0(k),y0(k),x0(el),y0(el));
         %dSqr = distVertexes(x(k),y(k),x(el),y(el));
         if (dSqr <= rangeRadiusSqr) == 1,
            % 根据节点通讯半径，划线
            xl = x(k):(x(el)-x(k))/10:x(el);
            yl = ((y(k) - y(el))/(x(k) - x(el)))*(xl - x(el)) + y(el);
            if get(handles.rbRed,'Value') == 1,
               plot(xl, yl, 'r');
            elseif get(handles.rbBlue,'Value') == 1,
               plot(xl, yl, 'b');
            elseif get(handles.rbWhite,'Value') == 1,
               plot(xl, yl, 'w');
            else
               plot(xl, yl, 'k');
            end    
         end
      end
   end
end


% --- Executes on button press in btnShowRadius.
function btnShowRadius_Callback(hObject, eventdata, handles)
% hObject    handle to btnShowRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global total;
global x;
global y;
global rangeRadius;

rangeRadius = str2double(get(handles.edtRangeRadius,'String'));

if get(handles.cbShowSingleRadius,'Value') == 1, 
   k = str2double(get(handles.edtVertexSeq,'String'));

   % 根据节点通讯半径，画圆
   alpha = 0:pi/10:2*pi; % 角度[0,2*pi] 
   xl = rangeRadius * cos(alpha) + x(k);
   yl = rangeRadius * sin(alpha) + y(k);
   if get(handles.rbRed,'Value') == 1,
      plot(xl, yl, 'r:');
   elseif get(handles.rbBlue,'Value') == 1,
      plot(xl, yl, 'b:');
   elseif get(handles.rbWhite,'Value') == 1,
      plot(xl, yl, 'w:');
   else
      plot(xl, yl, 'k:');
   end
else
   for k = 1:total,
      % 根据节点通讯半径，画圆
      alpha = 0:pi/10:2*pi; % 角度[0,2*pi] 
      xl = rangeRadius * cos(alpha) + x(k);
      yl = rangeRadius * sin(alpha) + y(k);
      if get(handles.rbRed,'Value') == 1,
         plot(xl, yl, 'r:');
      elseif get(handles.rbBlue,'Value') == 1,
         plot(xl, yl, 'b:');
      elseif get(handles.rbWhite,'Value') == 1,
         plot(xl, yl, 'w:');
      else
         plot(xl, yl, 'k:');
      end
   end
end


% --- Executes on button press in btnShowCutVertex.
function btnShowCutVertex_Callback(hObject, eventdata, handles)
% hObject    handle to btnShowCutVertex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global x;
global y;
global measMatrix;

cvSet = calcCutVertexSet(measMatrix);

if get(handles.rbRed,'Value') == 1,
   plot(x(cvSet), y(cvSet), 'r*');
elseif get(handles.rbBlue,'Value') == 1,
   plot(x(cvSet), y(cvSet), 'b*');
elseif get(handles.rbWhite,'Value') == 1,
   plot(x(cvSet), y(cvSet), 'w*');
else
   plot(x(cvSet), y(cvSet), 'k*');
end    


% --- Executes on button press in btnAdjMatrix.
function btnAdjMatrix_Callback(hObject, eventdata, handles)
% hObject    handle to btnAdjMatrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global total;
global x0;
global y0;
global x;
global y;
global rangeRadius;
global measMatrix;
global tdistMatrix;
global distMatrix;

disp('---------------------------------------------------------------------');
display(['Adjacent matrix begin time:', datestr(clock,0)]);

rangeRadius = str2double(get(handles.edtRangeRadius,'String'));
rangeRadiusSqr = rangeRadius^2;
tdistMatrix = zeros(total);
distMatrix = zeros(total);
measMatrix = zeros(total);

for i = 1:(total-1),
   for j = (i+1):total,
      dSqr = distVertexes(x0(i),y0(i),x0(j),y0(j));
      %dSqr = distVertexes(x(i),y(i),x(j),y(j));
      if (dSqr <= rangeRadiusSqr) == 1,
         %产生均值为0.6，方差为0.1的一个5*5的随机数方式
         %.6 + sqrt(0.1) * randn(5)
         tdistMatrix(i,j) = distVertexes(x0(i),y0(i),x0(j),y0(j));
         tdistMatrix(j,i) = tdistMatrix(i,j);
         distMatrix(i,j) = distVertexes(x(i),y(i),x(j),y(j));
         distMatrix(j,i) = distMatrix(i,j);
         measMatrix(i,j) = 1;
         measMatrix(j,i) = measMatrix(i,j);
      end
   end
end

save .\mat\distMatrix distMatrix;
save .\mat\measMatrix measMatrix;

display(['Adjacent matrix over time:', datestr(clock,0)]);
disp('---------------------------------------------------------------------');
   

% --- Executes on button press in btnDegreeSeq.
function btnDegreeSeq_Callback(hObject, eventdata, handles)
% hObject    handle to btnDegreeSeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%sort函数的调用格式:  
%sort(X)   功能：返回对向量X中的元素按列升序排列的新向量。
%[Y, I] = sort(A, dim, mode) 功能：对矩阵A的各列或各行重新排序，I记录Y中的元素在排序前A中位置，其中dim指明读A的列还是行进行排序。
%若dim=1，则按列排序；若dim=2，则按行排序。mode为排序的方式，取值'ascend'为升序，'descend'为降序。

global total;
global x;
global y;
global degreeSeq;
global vertexSeq;
global nbDegreeSeq;
global nbVertexSeq;
global measMatrix;

disp('---------------------------------------------------------------------');
display(['Degree seq begin time:', datestr(clock,0)]);

nbVertexDegree = zeros(total);
nbDegreeSeq = zeros(total);
nbVertexSeq = zeros(total);

vertexDeg = sum(measMatrix);
[degreeSeq, vertexSeq] = sort(vertexDeg, 2, 'descend');
display(['degree seq:', mat2str(degreeSeq)]);
display(['vertex seq:', mat2str(vertexSeq)]);
display(['average vertex degree:', num2str(mean(degreeSeq(:)))]);

for i = 1:total,
   for j = 1:total;
      nbVertexDegree(i,j) = measMatrix(i,j) * vertexDeg(j);
   end
  
   nbDegreeSeq(i,:) = sort(nbVertexDegree(i,:), 2, 'descend');
   nbVertexSeq(i,1:size(find(nbVertexDegree(i,:)>0),2)) = sort(find(nbVertexDegree(i,:)>0), 2, 'descend');
end

set(handles.tbSeq, 'data', [x; y; degreeSeq; vertexSeq]');

display(['Degree seq over time:', datestr(clock,0)]);
disp('---------------------------------------------------------------------');


% --- Executes on button press in btnTest.
function btnTest_Callback(hObject, eventdata, handles)
% hObject    handle to btnTest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global total;
global x0;
global y0;
global rangeRadius;
global rangeError;
global rangeShield;
global rangeVelocity;

display(['x0:', mat2str(x0)]);
display(['y0:', mat2str(y0)]);

rangeVelocity = 0;
rangeRadiusSqr = rangeRadius^2;

disp('---------------------------------------------------------------------');
display(['1 Adjacent matrix begin time:', datestr(clock,0)]);
         
% actual graph without error
measMatrix = zeros(total);
distLength = 0;
for i = 1:(total-1),
    for j = (i+1):total,
        dSqr = distVertexes(x0(i),y0(i),x0(j),y0(j));
        if (dSqr <= rangeRadiusSqr) == 1,
           distLength = distLength + 1;
                    
           measMatrix(i,j) = 1;
           measMatrix(j,i) = measMatrix(i,j);
                    
           xp = x0(i):(x0(j)-x0(i))/10:x0(j);
           yp = ((y0(i) - y0(j))/(x0(i) - x0(j)))*(xp - x0(j)) + y0(j);
           plot(xp, yp, 'k');
        end
    end
end
plot(x0, y0, 'k*');
        
display(['1 Adjacent matrix over time:', datestr(clock,0)]);
disp('---------------------------------------------------------------------');

disp('---------------------------------------------------------------------');
display(['2 Degree seq begin time:', datestr(clock,0)]);

nbVertexDegree = zeros(total);
nbDegreeSeq = zeros(total);
nbVertexSeq = zeros(total);

vertexDeg = sum(measMatrix);
[degreeSeq, vertexSeq] = sort(vertexDeg, 2, 'descend');
display(['degree seq:', mat2str(degreeSeq)]);
display(['vertex seq:', mat2str(vertexSeq)]);
display(['average vertex degree:', num2str(mean(degreeSeq(:)))]);

for i = 1:total,
    for j = 1:total;
        nbVertexDegree(i,j) = measMatrix(i,j) * vertexDeg(j);
    end
             
    nbDegreeSeq(i,:) = sort(nbVertexDegree(i,:), 2, 'descend');
    nbVertexSeq(i,1:size(find(nbVertexDegree(i,:)>0),2)) = sort(find(nbVertexDegree(i,:)>0), 2, 'descend');
end
         
display(['2 Degree seq over time:', datestr(clock,0)]);
disp('---------------------------------------------------------------------');

disp('---------------------------------------------------------------------');
display(['2.1 errorVector begin time:', datestr(clock,0)]);
         
%生成误差向量
rng('shuffle', 'v5normal');
errorVector = randn(1,distLength);
display(['errorVector:', mat2str(errorVector)]);
         
display(['2.1 errorVector over time:', datestr(clock,0)]);
disp('---------------------------------------------------------------------');
         
disp('---------------------------------------------------------------------');
display(['2.2 Bounded Distance begin time:', datestr(clock,0)]);

%调整误差，控制距离在Bounded内
tdistMatrix = zeros(total);
distMatrix = zeros(total);
idistLength = 0;
for i = 1:(total-1),
    for j = (i+1):total,
        dSqr = distVertexes(x0(i),y0(i),x0(j),y0(j));
        if (dSqr <= rangeRadiusSqr) == 1,
            idistLength = idistLength + 1;

            %产生均值为0.6，方差为0.1的一个5*5的随机数方式
            %.6 + sqrt(0.1) * randn(5)
            %tdistMatrix = true distance matrix
            tdistMatrix(i,j) = distVertexes(x0(i),y0(i),x0(j),y0(j));
            tdistMatrix(j,i) = tdistMatrix(i,j);

            %distMatrix = measurement distance matrix
            distMatrix(i,j) = (abs(sqrt(tdistMatrix(i,j)) + sqrt(rangeError)*errorVector(idistLength)))^2;
            if abs(sqrt(distMatrix(i,j)) - sqrt(tdistMatrix(i,j))) > rangeError*rangeShield,
                if sqrt(distMatrix(i,j)) > sqrt(tdistMatrix(i,j)),
                    distMatrix(i,j) = (sqrt(tdistMatrix(i,j))+rangeError*rangeShield)^2;
                elseif sqrt(distMatrix(i,j)) < sqrt(tdistMatrix(i,j)),
                    distMatrix(i,j) = (sqrt(tdistMatrix(i,j))-rangeError*rangeShield)^2;
                end
            end
            distMatrix(j,i) = distMatrix(i,j);
                    
            display(['i:', mat2str(i), ',j:', mat2str(j), ',distMatrix:', mat2str(sqrt(distMatrix(i,j))), ',tdistMatrix:', mat2str(sqrt(tdistMatrix(i,j))), ',diff:', mat2str(sqrt(distMatrix(i,j))-sqrt(tdistMatrix(i,j)))]);
        end
    end
end
            
display(['2.2 Bounded Distance over time:', datestr(clock,0)]);
disp('---------------------------------------------------------------------');

disp('---------------------------------------------------------------------');
display(['2.3 The initial triangle begin time:', datestr(clock,0)]);

%randomly initial triangle
triangle_n=0;
for i = 1:(total-2),
    for j = (i+1):(total-1),
        for k = (j+1):total,
            if (measMatrix(vertexSeq(i),vertexSeq(j)) == 1) && (measMatrix(vertexSeq(i),vertexSeq(k)) == 1) && (measMatrix(vertexSeq(j),vertexSeq(k)) == 1),
                 triangle_n = triangle_n + 1;
            end
        end
    end
end
            
triangle_r = unidrnd(triangle_n);
display(['triangle_r:', num2str(triangle_r)]);
            
triangle_i=0;
for i = 1:(total-2),
    for j = (i+1):(total-1),
        for k = (j+1):total,
            if (measMatrix(vertexSeq(i),vertexSeq(j)) == 1) && (measMatrix(vertexSeq(i),vertexSeq(k)) == 1) && (measMatrix(vertexSeq(j),vertexSeq(k)) == 1),
                triangle_i = triangle_i + 1;
                if triangle_i == triangle_r,
                    vertex_fo = vertexSeq(i); %o
                    vertex_fp = vertexSeq(j); %p
                    vertex_fq = vertexSeq(k); %q
         
                    display(['vertex_fo:', mat2str(vertex_fo), ' vertex_fp:', mat2str(vertex_fp), ' vertex_fq:', mat2str(vertex_fq)]);
                    plot(x0([vertex_fo,vertex_fp,vertex_fq]), y0([vertex_fo,vertex_fp,vertex_fq]), 'm*');

                    break;
                end
            end
        end
        if triangle_i == triangle_r,
            break;
        end
    end
    if triangle_i == triangle_r,
        break;
    end
end

distMatrix(vertex_fo,vertex_fp) = tdistMatrix(vertex_fo,vertex_fp);
distMatrix(vertex_fp,vertex_fo) = distMatrix(vertex_fo,vertex_fp);
distMatrix(vertex_fo,vertex_fq) = tdistMatrix(vertex_fo,vertex_fq);
distMatrix(vertex_fq,vertex_fo) = distMatrix(vertex_fo,vertex_fq);
distMatrix(vertex_fp,vertex_fq) = tdistMatrix(vertex_fp,vertex_fq);
distMatrix(vertex_fq,vertex_fp) = distMatrix(vertex_fp,vertex_fq);

display(['2.3 The first triangle over time:', datestr(clock,0)]);
disp('---------------------------------------------------------------------');

%AFALA
display(['rangeRadius:', num2str(rangeRadius), ' avgDegree:', num2str(mean(degreeSeq(:))), ' rangeError:', num2str(rangeError),' rangeShield:', num2str(rangeShield)]);
[xl, yl, vertexSeqLocl] = biTrilateration(0,0,total,x0,y0,vertex_fo,vertex_fp,vertex_fq,rangeRadius,rangeError,rangeVelocity,rangeShield,distMatrix,tdistMatrix,measMatrix,degreeSeq,vertexSeq,nbDegreeSeq,nbVertexSeq);


% --- Executes on button press in btnDeploy.
function btnDeploy_Callback(hObject, eventdata, handles)
% hObject    handle to btnDeploy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

networkset = [50,100,150,200];

axis([0 100 0 100]);
axis manual;
grid off; % grid on;
hold on; % hold off

for networkidx = 1:length(networkset);  
   for graphidx = 1:100,
      cla; %clear axes

      rng('shuffle', 'v5uniform');
      x0 = [unifrnd(0,100,networkset(networkidx),1)'];
      rng('shuffle', 'v5uniform');
      y0 = [unifrnd(0,100,networkset(networkidx),1)'];

      title('the graph of random node（x，y）');
      xlabel('X-axis');
      ylabel('Y-axis');
      plot(x0, y0, 'k*');

      save(['.\mat\x0_',num2str(networkidx),'_',num2str(graphidx),'.mat'],'x0');
      save(['.\mat\y0_',num2str(networkidx),'_',num2str(graphidx),'.mat'],'y0');
            
      hfigure = figure('visible', 'off'); %创建隐藏的窗口
      naxes = copyobj(handles.axes, hfigure); %将坐标轴区域复制到隐藏窗口
      set(naxes,'units','default','position','default');
      print(hfigure, '-dpng', ['.\mat\axes_',num2str(networkidx),'_',num2str(graphidx),'.png']); %输出到axespng图片
      delete(hfigure);
   end
end


% --- Executes on button press in btnAnalysis.
function btnAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to btnAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

format long;

networkset = [50,100,150,200];
radiusset = [20,30,40];
errorset =  [0,0.02,0.2,2];
shieldset = [1];
velocityset = [0:0.2:2];% represent velocity and direction as the distance change between nodes.

for networkidx = 1:length(networkset);  

   %for each graph
   for graphidx = 1:100;  

      cla; %clear axes
      clear global;
      axis([0 100 0 100]);
      axis manual;
      grid off; % grid on;
      hold on; % hold off
   
      total = networkset(networkidx);

      if exist(['./mat/x0_',num2str(networkidx),'_',num2str(graphidx),'.mat'],'file'),
         load(['./mat/x0_',num2str(networkidx),'_',num2str(graphidx),'.mat'],'x0');
      end
      if exist(['./mat/y0_',num2str(networkidx),'_',num2str(graphidx),'.mat'],'file'),
         load(['./mat/y0_',num2str(networkidx),'_',num2str(graphidx),'.mat'],'y0');
      end
   
      display(['x0:', mat2str(x0)]);
      display(['y0:', mat2str(y0)]);

      %for each radius
      for radiusidx = 1:length(radiusset);
      
         rangeRadius = radiusset(radiusidx);
         rangeRadiusSqr = radiusset(radiusidx)^2;

         disp('---------------------------------------------------------------------');
         display(['1 Adjacent matrix begin time:', datestr(clock,0)]);
         
%          figure(1);
%          cla; %clear axes
%          axis([0 100 0 100]);
%          axis manual;
%          grid off; % grid on;
%          hold on; % hold off

         % actual graph without error
         measMatrix = zeros(total);
         distLength = 0;
         for i = 1:(total-1),
            for j = (i+1):total,
               dSqr = distVertexes(x0(i),y0(i),x0(j),y0(j));
               if (dSqr <= rangeRadiusSqr) == 1,
                  distLength = distLength + 1;

                  measMatrix(i,j) = 1;
                  measMatrix(j,i) = measMatrix(i,j);
                    
                  xp = x0(i):(x0(j)-x0(i))/10:x0(j);
                  yp = ((y0(i) - y0(j))/(x0(i) - x0(j)))*(xp - x0(j)) + y0(j);
                  plot(xp, yp, 'k');
               end
            end
         end
         plot(x0, y0, 'k*');
        
         hfigure = figure('visible', 'off'); %创建隐藏的窗口
         naxes = copyobj(handles.axes, hfigure); %将坐标轴区域复制到隐藏窗口
         set(naxes,'units','default','position','default');
         print(hfigure, '-dpng', ['./mat/axes_',num2str(networkidx),'_',num2str(graphidx),'_',num2str(radiusidx),'_act.png']); %输出到axespng图片
         delete(hfigure);

         display(['1 Adjacent matrix over time:', datestr(clock,0)]);
         disp('---------------------------------------------------------------------');

         disp('---------------------------------------------------------------------');
         display(['2 Degree seq begin time:', datestr(clock,0)]);

         nbVertexDegree = zeros(total);
         nbDegreeSeq = zeros(total);
         nbVertexSeq = zeros(total);

         vertexDeg = sum(measMatrix);
         [degreeSeq, vertexSeq] = sort(vertexDeg, 2, 'descend');
         display(['degree seq:', mat2str(degreeSeq)]);
         display(['vertex seq:', mat2str(vertexSeq)]);
         display(['average vertex degree:', num2str(mean(degreeSeq(:)))]);

         for i = 1:total,
            for j = 1:total;
               nbVertexDegree(i,j) = measMatrix(i,j) * vertexDeg(j);
            end
             
            nbDegreeSeq(i,:) = sort(nbVertexDegree(i,:), 2, 'descend');
            nbVertexSeq(i,1:size(find(nbVertexDegree(i,:)>0),2)) = sort(find(nbVertexDegree(i,:)>0), 2, 'descend');
         end
         
         display(['2 Degree seq over time:', datestr(clock,0)]);
         disp('---------------------------------------------------------------------');
      
         %for each error
         for erroridx = 1:length(errorset);  %for each error
         
            rangeError = errorset(erroridx);
      
            disp('---------------------------------------------------------------------');
            display(['2.1 errorVector begin time:', datestr(clock,0)]);
         
            %生成误差向量
            rng('shuffle', 'v5normal');
            errorVector = randn(1,distLength);
            display(['errorVector:', mat2str(errorVector)]);
         
            display(['2.1 errorVector over time:', datestr(clock,0)]);
            disp('---------------------------------------------------------------------');
         
            %for each velocity
            for velocityidx = 1:length(velocityset);  %for each maxvelocity

               rangeVelocity = velocityset(velocityidx);
            
               disp('---------------------------------------------------------------------');
               display(['2.11 maxvelocityVector begin time:', datestr(clock,0)]);
         
               %生成速度向量
               rng('shuffle', 'v5uniform');
               velocityVector = rand(1,distLength);
               display(['velocityVector:', mat2str(velocityVector)]);
            
               display(['2.11 maxvelocityVector over time:', datestr(clock,0)]);
               disp('---------------------------------------------------------------------');
            
               %for each shield
               for shieldidx = 1:length(shieldset);  %for each shield
                  if shieldidx ~= 1,
                     continue;  
                  end
            
                  rangeShield = shieldset(shieldidx);
            
                  disp('---------------------------------------------------------------------');
                  display(['2.2 Bounded Distance begin time:', datestr(clock,0)]);

                  %调整误差，控制距离在Bounded内
                  tdistMatrix = zeros(total);
                  distMatrix = zeros(total);
                  idistLength = 0;
                  for i = 1:(total-1),
                     for j = (i+1):total,
                        dSqr = distVertexes(x0(i),y0(i),x0(j),y0(j));
                        if (dSqr <= rangeRadiusSqr) == 1,
                           idistLength = idistLength + 1;

                           %正态分布：产生均值为0.6，方差为0.1的一个5*5的随机数方式
                           %.6 + sqrt(0.1) * randn(5)
                           %均匀分布：产生(a,b) 区间上的均匀分布，只需对其做简单的线性变换即可：a+(b?a)?rand
                           %产生(?a,a)有可进一步化简为：?a+(a?(?a))?rand=a(2?rand?1)=(rand?1/2)?2?a
                           %(-5, 5)：-5+(5-(-5))*rand, (2*rand-1)*5
                           %tdistMatrix = true distance matrix
                           tdistMatrix(i,j) = distVertexes(x0(i),y0(i),x0(j),y0(j));
                           tdistMatrix(j,i) = tdistMatrix(i,j);

                           %distMatrix = measurement distance matrix
                           distMatrix(i,j) = (abs(sqrt(tdistMatrix(i,j)) + sqrt(rangeError)*errorVector(idistLength) + (velocityVector(idistLength) -0.5)*2*rangeVelocity))^2;
                           if abs(sqrt(distMatrix(i,j)) - sqrt(tdistMatrix(i,j))) > rangeError*rangeShield + rangeVelocity,
                              if sqrt(distMatrix(i,j)) > sqrt(tdistMatrix(i,j)),
                                 distMatrix(i,j) = (sqrt(tdistMatrix(i,j)) + rangeError*rangeShield + rangeVelocity)^2;
                              elseif sqrt(distMatrix(i,j)) < sqrt(tdistMatrix(i,j)),
                                 distMatrix(i,j) = (sqrt(tdistMatrix(i,j)) - rangeError*rangeShield - rangeVelocity)^2;
                              end
                           end
                           distMatrix(j,i) = distMatrix(i,j);
                    
                           display(['i:', mat2str(i), ',j:', mat2str(j), ',distMatrix:', mat2str(sqrt(distMatrix(i,j))), ',tdistMatrix:', mat2str(sqrt(tdistMatrix(i,j))), ',diff:', mat2str(sqrt(distMatrix(i,j))-sqrt(tdistMatrix(i,j)))]);
                        end
                     end
                  end
            
                  display(['2.2 Bounded Distance over time:', datestr(clock,0)]);
                  disp('---------------------------------------------------------------------');

                  disp('---------------------------------------------------------------------');
                  display(['2.3 The initial triangle begin time:', datestr(clock,0)]);

                  %randomly initial triangle
                  triangle_n=0;
                  for i = 1:(total-2),
                     for j = (i+1):(total-1),
                        for k = (j+1):total,
                           if (measMatrix(vertexSeq(i),vertexSeq(j)) == 1) && (measMatrix(vertexSeq(i),vertexSeq(k)) == 1) && (measMatrix(vertexSeq(j),vertexSeq(k)) == 1),
                              triangle_n = triangle_n + 1;
                           end
                        end
                     end
                  end
            
                  triangle_r = unidrnd(triangle_n);
                  display(['triangle_r:', num2str(triangle_r)]);
            
                  triangle_i=0;
                  for i = 1:(total-2),
                     for j = (i+1):(total-1),
                        for k = (j+1):total,
                           if (measMatrix(vertexSeq(i),vertexSeq(j)) == 1) && (measMatrix(vertexSeq(i),vertexSeq(k)) == 1) && (measMatrix(vertexSeq(j),vertexSeq(k)) == 1),
                              triangle_i = triangle_i + 1;
                              if triangle_i == triangle_r,
                                 vertex_fo = vertexSeq(i); %o
                                 vertex_fp = vertexSeq(j); %p
                                 vertex_fq = vertexSeq(k); %q
         
                                 display(['vertex_fo:', mat2str(vertex_fo), ' vertex_fp:', mat2str(vertex_fp), ' vertex_fq:', mat2str(vertex_fq)]);
                                 plot(x0([vertex_fo,vertex_fp,vertex_fq]), y0([vertex_fo,vertex_fp,vertex_fq]), 'm*');

                                 break;
                              end
                           end
                        end
                        if triangle_i == triangle_r,
                           break;
                        end
                     end
                     if triangle_i == triangle_r,
                        break;
                     end
                  end
            
%                   %speciallty initial triangle
%                   flag=0;
%                   for i = 1:(total-2),
%                      for j = (i+1):(total-1),
%                         for k = (j+1):total,
%                            if (measMatrix(vertexSeq(i),vertexSeq(j)) == 1) && (measMatrix(vertexSeq(i),vertexSeq(k)) == 1) && (measMatrix(vertexSeq(j),vertexSeq(k)) == 1)...
%                               && sqrt(distMatrix(vertexSeq(i),vertexSeq(j))) > rangeRadius*4/5 && sqrt(distMatrix(vertexSeq(i),vertexSeq(k))) > rangeRadius*4/5 && sqrt(distMatrix(vertexSeq(j),vertexSeq(k))) > rangeRadius*4/5 ...
%                               && sqrt(distMatrix(vertexSeq(i),vertexSeq(j))) > 2*rangeError*rangeShield && sqrt(distMatrix(vertexSeq(i),vertexSeq(k))) > 2*rangeError*rangeShield && sqrt(distMatrix(vertexSeq(j),vertexSeq(k))) > 2*rangeError*rangeShield,
%                                vertex_fo = vertexSeq(i); %o
%                                vertex_fp = vertexSeq(j); %p
%                                vertex_fq = vertexSeq(k); %q
%          
%                                display(['vertex_fo:', mat2str(vertex_fo), ' vertex_fp:', mat2str(vertex_fp), ' vertex_fq:', mat2str(vertex_fq)]);
%                                plot(x0([vertex_fo,vertex_fp,vertex_fq]), y0([vertex_fo,vertex_fp,vertex_fq]), 'm*');
% 
%                                flag=1;
%                                break;
%                            end
%                         end
%                         if flag == 1,
%                            break;
%                         end
%                      end
%                      if flag == 1,
%                         break;
%                      end
%                   end

                  distMatrix(vertex_fo,vertex_fp) = tdistMatrix(vertex_fo,vertex_fp);
                  distMatrix(vertex_fp,vertex_fo) = distMatrix(vertex_fo,vertex_fp);
                  distMatrix(vertex_fo,vertex_fq) = tdistMatrix(vertex_fo,vertex_fq);
                  distMatrix(vertex_fq,vertex_fo) = distMatrix(vertex_fo,vertex_fq);
                  distMatrix(vertex_fp,vertex_fq) = tdistMatrix(vertex_fp,vertex_fq);
                  distMatrix(vertex_fq,vertex_fp) = distMatrix(vertex_fp,vertex_fq);

                  display(['2.3 The first triangle over time:', datestr(clock,0)]);
                  disp('---------------------------------------------------------------------');
         
                  display(['rangeRadius:', num2str(rangeRadius), ' avgDegree:', num2str(mean(degreeSeq(:))), ' rangeError:', num2str(rangeError), ' rangeVelocity:', num2str(rangeVelocity),' rangeShield:', num2str(rangeShield)]);
                  [xl, yl, vertexSeqLocl] = biTrilateration(networkidx,graphidx,total,x0,y0,vertex_fo,vertex_fp,vertex_fq,rangeRadius,rangeError,rangeVelocity,rangeShield,distMatrix,tdistMatrix,measMatrix,degreeSeq,vertexSeq,nbDegreeSeq,nbVertexSeq);

%                   hfigure = figure('visible', 'off'); %创建隐藏的窗口
%                   naxes = copyobj(handles.axes, hfigure); %将坐标轴区域复制到隐藏窗口
%                   set(naxes,'units','default','position','default');
%                   print(hfigure, '-dpng', ['./mat/axes_',num2str(networkidx),'_',num2str(graphidx),'_',num2str(radiusidx),'_',num2str(erroridx),'_',num2str(velocityidx),'_',num2str(shieldidx),'_bitri_1.png']); %输出到axespng图片
%                   delete(hfigure);

                  if rangeError == 0,
                     break;
                  end
               end
            end
         end
      end
   end
end


% --- Executes on button press in btnSavePng.
function btnSavePng_Callback(hObject, eventdata, handles)
% hObject    handle to btnSavePng (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hfigure = figure('visible', 'off'); %创建隐藏的窗口
naxes = copyobj(handles.axes, hfigure); %将坐标轴区域复制到隐藏窗口
set(naxes,'units','default','position','default');
print(hfigure, '-dpng', '.\mat\axes.png'); %输出到axespng图片
delete(hfigure);

   
% --- Executes on button press in btnClose.
function btnClose_Callback(hObject, eventdata, handles)
% hObject    handle to btnClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close;
