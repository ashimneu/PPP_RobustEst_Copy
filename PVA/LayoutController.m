% Notes: This script was developed to arrange the figure windows generated 
% by the function plotTraj on multiple screens. 
% There are two different layouts I have created for the figure window .

% Author: Ashim Neupane; Email: aneup001@ucr.edu; Last Updated: 4/1/2021

solvername = "kf"; % choices: "kf", "lts", "td", "raps","mshb","mstk"
option.eb_LTS  = output.p.eb_LTS; 
option.eb_RAPS = output.p.eb_RAPS; 
option.eb_TD   = output.p.eb_TD;
option.eb_MShb = output.p.eb_MShb;
option.eb_MStk = output.p.eb_MStk;

loadgui(solvername,option)

function loadgui(solvername,option)
% For KF Plots
solver_Plots = createFigArray(solvername);

Ctrlboxtag = strcat("FigLayoutControlBox");
fignum = getfignum(Ctrlboxtag);    
figbox = figure(fignum); clf; hold on
set(figbox,'Name','Solver Figures Layout Controller','NumberTitle','off');

% Create different panel to hold the buttons & text boxes
CtrlPanel = uipanel(figbox, 'Position', [0 0 1 1], 'BackgroundColor', [1 1 1]);

btnWidth  = 0.15;
btnHeight = 0.05;

%%%%%%%%%%%%%%%%%%%%%% Stacked Layout (Layout - 1) %%%%%%%%%%%%%%%%%%%%%%%%
lay1.Title_txt = uicontrol(CtrlPanel, 'Style','text', 'String','Stacked Layout','Units', 'normalized', 'Position', [0.05 0.90 0.2 0.05]);

% Screen Text & Textboxes
lay1.Screen_A_y = 0.85; lay1.Screen_A_x = 0.05;
lay1.txt_A = uicontrol(CtrlPanel,'Style','text', 'String','Screen A',...
    'Units', 'normalized', 'Position', [lay1.Screen_A_x lay1.Screen_A_y btnWidth btnHeight]);
lay1.txt_fwidth_A = uicontrol(CtrlPanel,'Style','text', 'String','Fig Width',...
    'Units', 'normalized', 'Position', [(lay1.Screen_A_x+btnWidth) lay1.Screen_A_y btnWidth btnHeight]);
lay1.etxt_fwidth_A = uicontrol('style','Edit','string','-','Units','normalized',...
                   'Position',[(lay1.Screen_A_x+2*btnWidth) lay1.Screen_A_y btnWidth btnHeight],'backgroundcolor','w','Tag','EditField');
lay1.txt__fOdx_A = uicontrol(CtrlPanel, 'Style','text', 'String','Left edge dist.',...
    'Units', 'normalized', 'Position', [(lay1.Screen_A_x+3*btnWidth) lay1.Screen_A_y btnWidth btnHeight]);               
lay1.etxt_fOdx_A = uicontrol('style','Edit','string','-','Units','normalized',...
   'Position',[(lay1.Screen_A_x+4*btnWidth) lay1.Screen_A_y btnWidth btnHeight],'backgroundcolor','w','Tag','EditField');

lay1.Btny = 0.80;
lay1.BtnHeight = 0.05;
lay1.Scn1Btn = uicontrol(CtrlPanel, 'Style', 'pushbutton', ...
    'String', 'Screen 1', 'Units', 'normalized', 'Position', [(lay1.Screen_A_x+0*btnWidth) lay1.Btny btnWidth lay1.BtnHeight], ...
    'Callback', @(src, event) lay1BtnPressCbak(solver_Plots,1,lay1));
lay1.Scn2Btn = uicontrol(CtrlPanel, 'Style', 'pushbutton', ...
    'String', 'Screen 2', 'Units', 'normalized', 'Position', [(lay1.Screen_A_x+1*btnWidth) lay1.Btny btnWidth lay1.BtnHeight], ...
    'Callback', @(src, event) lay1BtnPressCbak(solver_Plots,2,lay1));
lay1.Scn3Btn = uicontrol(CtrlPanel, 'Style', 'pushbutton', ...
    'String', 'Screen 3', 'Units', 'normalized', 'Position', [(lay1.Screen_A_x+2*btnWidth) lay1.Btny btnWidth lay1.BtnHeight], ...
    'Callback', @(src, event) lay1BtnPressCbak(solver_Plots,3,lay1));
% set(KF_Plots(2),'units','normalized','position',[1.25 -0.20 0.280 1.1])
lay1.ScnBtnUp = uicontrol(CtrlPanel, 'Style', 'pushbutton', ...
    'String', 'Up', 'Units', 'normalized', 'Position', [(lay1.Screen_A_x+3*btnWidth) lay1.Btny btnWidth lay1.BtnHeight], ...
    'Callback', @(src, event) WindowStateBtnPressCbak(solver_Plots,"up"));
lay1.ScnBtnDown = uicontrol(CtrlPanel, 'Style', 'pushbutton', ...
    'String', 'Down', 'Units', 'normalized', 'Position', [(lay1.Screen_A_x+4*btnWidth) lay1.Btny btnWidth lay1.BtnHeight], ...
    'Callback', @(src, event) WindowStateBtnPressCbak(solver_Plots,"down"));

%%%%%%%%%%%%%%%%% Split Screen Layout (Layout - 2) %%%%%%%%%%%%%%%%%%%%%%%%
lay2.Title_txt = uicontrol(CtrlPanel, 'Style','text', 'String','Two Screen Layout','Units', 'normalized', 'Position', [0.05 0.60 0.2 0.05]);
lay2.Screen_A_x = 0.05;

% Screen A Text & Textboxes
lay2.Screen_A_y = 0.55;
lay2.txt_A = uicontrol(CtrlPanel,'Style','text', 'String','Screen A',...
    'Units', 'normalized', 'Position', [lay2.Screen_A_x lay2.Screen_A_y btnWidth btnHeight]);
lay2.txt_fwidth_A = uicontrol(CtrlPanel,'Style','text', 'String','Fig Width',...
    'Units', 'normalized', 'Position', [(lay2.Screen_A_x+btnWidth) lay2.Screen_A_y btnWidth btnHeight]);
lay2.etxt_fwidth_A = uicontrol('style','Edit','string','-','Units','normalized',...
                   'Position',[(lay2.Screen_A_x+2*btnWidth) lay2.Screen_A_y btnWidth btnHeight],'backgroundcolor','w','Tag','EditField');
lay2.txt__fOdx_A = uicontrol(CtrlPanel, 'Style','text', 'String','Left edge dist.',...
    'Units', 'normalized', 'Position', [(lay2.Screen_A_x+3*btnWidth) lay2.Screen_A_y btnWidth btnHeight]);               
lay2.etxt_fOdx_A = uicontrol('style','Edit','string','-','Units','normalized',...
   'Position',[(lay2.Screen_A_x+4*btnWidth) lay2.Screen_A_y btnWidth btnHeight],'backgroundcolor','w','Tag','EditField');

% Screen B Text & Textboxes
lay2.Screen_B_y = 0.5;
lay2.txt_Screen_B = uicontrol(CtrlPanel,'Style','text', 'String','Screen B',...
    'Units', 'normalized', 'Position', [lay2.Screen_A_x lay2.Screen_B_y btnWidth btnHeight]);
lay2.txt_fwidth_B = uicontrol(CtrlPanel,'Style','text', 'String','Fig Width',...
    'Units', 'normalized', 'Position', [(lay2.Screen_A_x+btnWidth) lay2.Screen_B_y btnWidth btnHeight]);
lay2.etxt_fwidth_B = uicontrol('style','Edit','string','-','Units','normalized',...
                   'Position',[(lay2.Screen_A_x+2*btnWidth) lay2.Screen_B_y btnWidth btnHeight],'backgroundcolor','w','Tag','EditField');
lay2.txt_fOdx_B = uicontrol(CtrlPanel, 'Style','text', 'String','Left edge dist.',...
    'Units', 'normalized', 'Position', [(lay2.Screen_A_x+3*btnWidth) lay2.Screen_B_y btnWidth btnHeight]);               
lay2.etxt_fOdx_B = uicontrol('style','Edit','string','-','Units','normalized',...
   'Position',[(lay2.Screen_A_x+4*btnWidth) lay2.Screen_B_y btnWidth btnHeight],'backgroundcolor','w','Tag','EditField');

% Buttons
lay2_Btny = 0.45;
lay2.Scn12Btn = uicontrol(CtrlPanel, 'Style', 'pushbutton', ...
    'String', 'Screen 1,2', 'Units', 'normalized', 'Position', [lay2.Screen_A_x lay2_Btny btnWidth btnHeight], ...
    'Callback', @(src, event) lay2BtnPressCbak(solver_Plots,1,2,lay2));
lay2.Scn23Btn = uicontrol(CtrlPanel, 'Style', 'pushbutton', ...
    'String', 'Screen 2,3', 'Units', 'normalized', 'Position', [(lay2.Screen_A_x+btnWidth) lay2_Btny btnWidth btnHeight], ...
    'Callback', @(src, event) lay2BtnPressCbak(solver_Plots,2,3,lay2));
lay2.Scn13Btn = uicontrol(CtrlPanel, 'Style', 'pushbutton', ...
    'String', 'Screen 1,3', 'Units', 'normalized', 'Position', [(lay2.Screen_A_x+2*btnWidth) lay2_Btny btnWidth btnHeight], ...
    'Callback', @(src, event) lay2BtnPressCbak(solver_Plots,1,3,lay2));
lay2.ScnABtnUp = uicontrol(CtrlPanel, 'Style', 'pushbutton', ...
    'String', 'Up', 'Units', 'normalized', 'Position', [(lay2.Screen_A_x+3*btnWidth) lay2_Btny btnWidth btnHeight], ...
    'Callback', @(src, event) WindowStateBtnPressCbak(solver_Plots,"up"));
lay2.ScnABtnDown = uicontrol(CtrlPanel, 'Style', 'pushbutton', ...
    'String', 'Down', 'Units', 'normalized', 'Position', [(lay2.Screen_A_x+4*btnWidth) lay2_Btny btnWidth btnHeight], ...
    'Callback', @(src, event) WindowStateBtnPressCbak(solver_Plots,"down"));

lay2_Btny = 0.40;
lay2.Scn21Btn = uicontrol(CtrlPanel, 'Style', 'pushbutton', ...
    'String', 'Screen 2,1', 'Units', 'normalized', 'Position', [lay2.Screen_A_x lay2_Btny btnWidth btnHeight], ...
    'Callback', @(src, event) lay2BtnPressCbak(solver_Plots,2,1,lay2));
lay2.Scn32Btn = uicontrol(CtrlPanel, 'Style', 'pushbutton', ...
    'String', 'Screen 3,2', 'Units', 'normalized', 'Position', [(lay2.Screen_A_x+btnWidth) lay2_Btny btnWidth btnHeight], ...
    'Callback', @(src, event) lay2BtnPressCbak(solver_Plots,3,2,lay2));
lay2.Scn31Btn = uicontrol(CtrlPanel, 'Style', 'pushbutton', ...
    'String', 'Screen 3,1', 'Units', 'normalized', 'Position', [(lay2.Screen_A_x+2*btnWidth) lay2_Btny btnWidth btnHeight], ...
    'Callback', @(src, event) lay2BtnPressCbak(solver_Plots,3,1,lay2));
lay2.ScnBBtnUp = uicontrol(CtrlPanel, 'Style', 'pushbutton', ...
    'String', 'Up', 'Units', 'normalized', 'Position', [(lay2.Screen_A_x+3*btnWidth) lay2_Btny btnWidth btnHeight], ...
    'Callback', @(src, event) WindowStateBtnPressCbak(solver_Plots,"up"));
lay2.ScnBBtnDown = uicontrol(CtrlPanel, 'Style', 'pushbutton', ...
    'String', 'Down', 'Units', 'normalized', 'Position', [(lay2.Screen_A_x+4*btnWidth) lay2_Btny btnWidth btnHeight], ...
    'Callback', @(src, event) WindowStateBtnPressCbak(solver_Plots,"down"));

%%%%%%%%%%%%%%%%%%%%%%%%%% Bottom Section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BtmSec_y = 0.3;
BtmSec.Title_txt = uicontrol(CtrlPanel, 'Style','text', 'String','Algorithm Selected','Units', 'normalized', 'Position', [0.05 BtmSec_y 0.2 0.05]);
txtWidth = 0.12;

% Screen Text & Textboxes
BtmSec.algo_txt_y = 0.25; BtmSec.algo_txt_x = 0.05;
BtmSec.Btn_KF = uicontrol(CtrlPanel,'Style','pushbutton', 'String','KF',...
    'Units', 'normalized', 'Position', [BtmSec.algo_txt_x BtmSec.algo_txt_y txtWidth btnHeight]);
BtmSec.Btn_LTS = uicontrol(CtrlPanel,'Style','pushbutton', 'String','LTS',...
    'Units', 'normalized', 'Position', [(BtmSec.algo_txt_x + txtWidth) BtmSec.algo_txt_y txtWidth btnHeight]);
BtmSec.Btn_TD = uicontrol(CtrlPanel,'Style','pushbutton', 'String','TD',...
    'Units', 'normalized', 'Position', [(BtmSec.algo_txt_x + 2*txtWidth) BtmSec.algo_txt_y txtWidth btnHeight]);
BtmSec.Btn_RAPS = uicontrol(CtrlPanel,'Style','pushbutton', 'String','RAPS',...
    'Units', 'normalized', 'Position', [(BtmSec.algo_txt_x + 3*txtWidth) BtmSec.algo_txt_y txtWidth btnHeight]);
BtmSec.Btn_MShb = uicontrol(CtrlPanel,'Style','pushbutton', 'String','MShb',...
    'Units', 'normalized', 'Position', [(BtmSec.algo_txt_x + 4*txtWidth) BtmSec.algo_txt_y txtWidth btnHeight]);
BtmSec.Btn_MStk = uicontrol(CtrlPanel,'Style','pushbutton', 'String','MStk',...
    'Units', 'normalized', 'Position', [(BtmSec.algo_txt_x + 5*txtWidth) BtmSec.algo_txt_y txtWidth btnHeight]);
set(BtmSec.Btn_KF,'Callback', @(src, event) SelectionBtnPressCbak("kf",BtmSec,option));
set(BtmSec.Btn_LTS,'Callback', @(src, event) SelectionBtnPressCbak("lts",BtmSec,option));
set(BtmSec.Btn_RAPS,'Callback', @(src, event) SelectionBtnPressCbak("raps",BtmSec,option));
set(BtmSec.Btn_TD,'Callback', @(src, event) SelectionBtnPressCbak("td",BtmSec,option));
set(BtmSec.Btn_MShb,'Callback', @(src, event) SelectionBtnPressCbak("mshb",BtmSec,option));
set(BtmSec.Btn_MStk,'Callback', @(src, event) SelectionBtnPressCbak("mstk",BtmSec,option));

if solvername == "kf"
    set(BtmSec.Btn_KF,'BackgroundColor','g');
        
elseif solvername == "lts" & option.eb_LTS == 1        
    set(BtmSec.Btn_LTS,'BackgroundColor','g');
            
elseif solvername == "raps" & option.eb_RAPS == 1        
    set(BtmSec.Btn_RAPS,'BackgroundColor','g');
           
elseif solvername == "td" & option.eb_TD == 1         %#ok<*AND2>
    set(BtmSec.Btn_TD,'BackgroundColor','g');

elseif solvername == "mshb" & option.eb_MShb == 1        
    set(BtmSec.Btn_MShb,'BackgroundColor','g');

elseif solvername == "mstk" & option.eb_MStk == 1        
    set(BtmSec.Btn_MStk,'BackgroundColor','g');
end

end
%--------------------------------------------------------------------------
function WindowStateBtnPressCbak(fig_Plots,desired_state)
    fig_Plots = removeClosedFigs(fig_Plots);    
    if lower(desired_state) == "down"
        set(fig_Plots,'WindowState','minimized');
    elseif lower(desired_state) == "up"
        set(fig_Plots,'WindowState','normal');
    end
 
end
%--------------------------------------------------------------------------
function lay2BtnPressCbak(fig_Plots,screen_A,screen_B,lay2)
    % Call Back fuction that arranges the figure windows as defined in
    % Layout 2.
    fig_Plots = removeClosedFigs(fig_Plots);
    
    fig_width_A = str2double(lay2.etxt_fwidth_A.String);
    fig_width_B = str2double(lay2.etxt_fwidth_B.String);
    fig_Origin_dx_A = str2double(lay2.etxt_fOdx_A.String);
    fig_Origin_dx_B = str2double(lay2.etxt_fOdx_B.String);
    if isnan(fig_width_A), fig_width_A  = getParam(screen_A,2).fig_width; insertVal(lay2.etxt_fwidth_A,fig_width_A); end
    if isnan(fig_width_B), fig_width_B  = getParam(screen_B,2).fig_width; insertVal(lay2.etxt_fwidth_B,fig_width_A); end
    if isnan(fig_Origin_dx_A), fig_Origin_dx_A = getParam(screen_A,2).fig_Origin_dx; insertVal(lay2.etxt_fOdx_A,fig_Origin_dx_A); end
    if isnan(fig_Origin_dx_B), fig_Origin_dx_B = getParam(screen_B,2).fig_Origin_dx; insertVal(lay2.etxt_fOdx_B,fig_Origin_dx_B); 
    end

    fig_Origin_x_A  = 0;    
    fig_Origin_x_B  = 0;    
    
    for i = 1:numel(fig_Plots)
        if i<=5
            if ~isempty(fig_Plots(i))
                set(fig_Plots(i),'units','normalized','outerposition',ZerothFigPosition(screen_A,fig_width_A)+ [fig_Origin_x_A 0 0 0])
                fig_Origin_x_A = fig_Origin_x_A + fig_Origin_dx_A;
            end            
        else
            if ~isempty(fig_Plots(i))
                set(fig_Plots(i),'units','normalized','outerposition',ZerothFigPosition(screen_B,fig_width_B)+ [fig_Origin_x_B 0 0 0])
                fig_Origin_x_B = fig_Origin_x_B + fig_Origin_dx_B;
            end            
        end
    end

end
function insertVal(handle,value)
    handle.String = value;
end
function param = getParam(screen,layout)
% Default Parameters
    if layout== 1
        switch screen
            case 1
                param.fig_Origin_x = 0;    
                param.fig_Origin_y = 0.025;
                param.fig_width    = 0.3; % default width  
                param.fig_height   = 1 - param.fig_Origin_y;
                param.fig_Origin_dx = 0.1;
            case 2
                param.fig_Origin_x = 1.25;    
                param.fig_Origin_y = -0.2;
                param.fig_width   = 0.3; % default width   
                param.fig_height  = 1.2;
                param.fig_Origin_dx = 0.15;
            case 3
                param.fig_Origin_x = 2.5;     
                param.fig_Origin_y = -0.2;
                param.fig_width   = 0.3; % default width   
                param.fig_height  = 1.2;
                param.fig_Origin_dx = 0.15;
        end
    elseif layout == 2
        switch screen
            case 1
                param.fig_Origin_x = 0;    
                param.fig_Origin_y = 0.025;
                param.fig_width    = 0.3; % default width  
                param.fig_height   = 1 - param.fig_Origin_y;
                param.fig_Origin_dx = 0.175;
            case 2
                param.fig_Origin_x = 1.25;    
                param.fig_Origin_y = -0.2;
                param.fig_width   = 0.3; % default width   
                param.fig_height  = 1.2;
                param.fig_Origin_dx = 0.25;
            case 3
                param.fig_Origin_x = 2.5;     
                param.fig_Origin_y = -0.2;
                param.fig_width   = 0.3; % default width   
                param.fig_height  = 1.2;
                param.fig_Origin_dx = 0.25;
        end
    end
end
function figPosition = ZerothFigPosition(screen,fig_width_in)
    % 0th Place for a figure in a screen is defined as the left most side of the chosen screen. Here default params are used for this purpose.
    param = getParam(screen,2);
    if ~isempty(fig_width_in)
        param.fig_width = fig_width_in;
    end
    figPosition = [param.fig_Origin_x param.fig_Origin_y param.fig_width param.fig_height];
end
function lay1BtnPressCbak(fig_Plots,screen_A,lay1)    
    % Call Back fuction that arranges the figure windows in screen 1. Layout 1 is used.
    fig_Plots = removeClosedFigs(fig_Plots);
    
    fig_width_A = str2double(lay1.etxt_fwidth_A.String);
    fig_Origin_dx_A = str2double(lay1.etxt_fOdx_A.String);
    if isnan(fig_width_A), fig_width_A = getParam(screen_A,1).fig_width; insertVal(lay1.etxt_fwidth_A,fig_width_A); end
    if isnan(fig_Origin_dx_A), fig_Origin_dx_A = getParam(screen_A,1).fig_Origin_dx; insertVal(lay1.etxt_fOdx_A,fig_Origin_dx_A); end
        
%     switch screen_A
%         case 1
%             Originx  = 0;       Originy   = 0.025;
%             figwidth = 0.3;     figheight = 1 - Originy;
%             dOriginx = 0.1;     
%         case 2
%             Originx  = 1.25;    Originy   = -0.2;
%             figwidth = 0.28;    figheight = 1.2;
%             dOriginx = 0.1;
%         case 3
%             Originx  = 2.5;     Originy   = -0.2;
%             figwidth = 0.28;    figheight = 1.2;
%             dOriginx = 0.1;
%     end
%         dOriginx = 0.12;
%         figwidth = 0.28;
        
        
%     for i = 1:numel(fig_Plots)
%         if ~isempty(fig_Plots(i))
%             set(fig_Plots(i),'units','normalized','outerposition',[Originx Originy figwidth figheight])
%         end
%         Originx = Originx + dOriginx;
%     end
    fig_Origin_x_A  = 0; 
    for i = 1:numel(fig_Plots)
            if ~isempty(fig_Plots(i))
                set(fig_Plots(i),'units','normalized','outerposition',ZerothFigPosition(screen_A,fig_width_A)+ [fig_Origin_x_A 0 0 0])
                fig_Origin_x_A = fig_Origin_x_A + fig_Origin_dx_A;
            end        
    end

end
%--------------------------------------------------------------------------
function SelectionBtnPressCbak(solvername,BtmSec,option)

BtnArray = [BtmSec.Btn_KF;BtmSec.Btn_LTS;BtmSec.Btn_RAPS;BtmSec.Btn_TD];

if solvername == "kf"
    setDefaultBtnColor(BtnArray);
    set(BtmSec.Btn_KF,'BackgroundColor','g');
    solver_Plots = createFigArray(solvername);
    set(BtmSec.Btn_KF,'UserData',solver_Plots);
    
elseif solvername == "lts" & option.eb_LTS == 1 
    setDefaultBtnColor(BtnArray);
    set(BtmSec.Btn_LTS,'BackgroundColor','g');
    solver_Plots = createFigArray(solvername);
    set(BtmSec.Btn_KF,'UserData',solver_Plots);
        
elseif solvername == "raps" & option.eb_RAPS == 1
    setDefaultBtnColor(BtnArray);
    set(BtmSec.Btn_RAPS,'BackgroundColor','g');
    solver_Plots = createFigArray(solvername);
    set(BtmSec.Btn_KF,'UserData',solver_Plots);
        
elseif solvername == "td" & option.eb_TD == 1         %#ok<*AND2>
    setDefaultBtnColor(BtnArray);
    set(BtmSec.Btn_TD,'BackgroundColor','g');
    solver_Plots = createFigArray(solvername);
    set(BtmSec.Btn_KF,'UserData',solver_Plots);        
else
%     fprintf('Warning! No figure for the algorithm you selected was found. \n')
    errordlg('Warning! No figures for the mesurement selection algorithm you selected were found.','Figures not found.')
    return
end
    loadgui(solvername,option)

end
function setDefaultBtnColor(BtnArray)
    BtnArray = BtnArray(isvalid(BtnArray));
    set(BtnArray,'BackgroundColor',[0.94 0.94 0.94]); % set to default grayish color
end
%--------------------------------------------------------------------------
function fig_array = createFigArray(solvername)
    fig_pos = findobj('type','figure','tag',strcat("pos_",solvername));
    fig_vel = findobj('type','figure','tag',strcat("vel_",solvername));
    fig_acl = findobj('type','figure','tag',strcat("acl_",solvername));
    fig_clk = findobj('type','figure','tag',strcat("clk_biasdrift_",solvername));
    fig_isb = findobj('type','figure','tag',strcat("ISB_",solvername));
    fig_res = findobj('type','figure','tag',strcat("res_",solvername));
    fig_max_res = findobj('type','figure','tag',strcat("max_res_",solvername));
    fig_infopos = findobj('type','figure','tag',strcat("info_pos_",solvername));
    fig_infoclk = findobj('type','figure','tag',strcat("info_clk_",solvername));
    fig_pos_err = findobj('type','figure','tag',strcat("pos_err_",solvername));    
    fig_estpred_cov = findobj('type','figure','tag',strcat("estpred_cov_",solvername));
    
    fig_array = [fig_pos; fig_vel; fig_acl; fig_clk; fig_isb; fig_res;...
        fig_max_res; fig_infopos; fig_infoclk; fig_pos_err; fig_estpred_cov];
end
%--------------------------------------------------------------------------
function fig_array_out = removeClosedFigs(fig_array)
    fig_array_out = fig_array(isvalid(fig_array));
end