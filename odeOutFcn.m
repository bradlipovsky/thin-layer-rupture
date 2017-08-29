function stop = odeOutFcn( t, ~, state , varargin )
% Stop intergration based on some criterion.
stop = false;
persistent step
persistent t0
persistent tt0
spacing = 1e3;


switch state
    case []
        % several output values may be calculated in one major time
        % step
        step = step + numel(t); 

%         if step > varargin{1}.MaxTimeSteps;
%             stop = true;
%         end
        if mod(step,spacing)==0
            tt=toc;
            dt = (t-t0)/spacing;
            dtt = (tt-tt0);
            disp(['Step ' num2str(step) ', wall ' num2str(tt) ', sim ' num2str(t)...
            	  ', s/w ' num2str(dt/dtt) ', avg step ' num2str(dt)]);
            t0 = t;
            tt0 = tt;
        end

    case 'init'
        step = 0;
        t0=0;
        tt0=0;
        disp(' ');
end
end