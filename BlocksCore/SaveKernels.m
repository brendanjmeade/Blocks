function SaveKernels(Partials, Command, varargin)


if strcmp(Command.saveKernels, 'yes') == 1;
   if exist('tempkernels.mat', 'file')
      for i = 1:length(varargin)
         save('tempkernels.mat', '-struct', 'Partials', varargin{i}, '-append');
      end
   else
      save('tempkernels.mat', '-struct', 'Partials');
   end
end