function [mmag, tarea] = couprelease(p, coup, area, incr, wcscale, subset)
% couprelease   Earthquake magnitude releasing coupling fraction.
%   couprelease(P, COUP, AREA, INCR) calculates the earthquake magnitude
%   necessary to release coupling accumulated within the areas with
%   a specified range of coupling fractions. P is a patch structure
%   containing the network of triangular elements on which the 
%   coupling distribution has been estimated. COUP is an n-by-1 array
%   giving the actual coupling distribution, AREA is an n-by-1 array
%   giving the areas of the elements on which the coupling distribution
%   is estimated, and INCR is a vector giving the coupling fractions
%   for which the required earthquake magnitude should be calculated.
%   Magnitude is calculated using the Wells and Coppersmith (1994)
%   scaling relations between rupture area and moment magnitude.
%
%   couprepease(COUP, AREA, INCR, COEFF) permits specification of the
%   Wells and Coppersmith scaling coefficients as the 4-element vector
%   COEFF, with stucture COEFF = [A B Auncertainty Buncertainty].
%
%   MMAG = couprelease(...) returns the moment magnitudes to MMAG.
%
%   [MMAG, TAREA] = couprelease(...) also returns the total areas in 
%   square kilometers.
%
%   The idea is that, for a given value of INCR, all elements of that
%   coupling fraction or greater are found and their areas are summed,
%   and returned as TAREA. TAREA is also converted to moment magnitude
%   using the Wells and Coppersmith "All" relationship of Table 2A.
%
%

% Check increment difference
id = mean(diff(incr));

% Allocate space for total area
tarea = NaN(sum(p.nEl), length(incr));

% For each coupling fraction,
for i = 1:length(incr)
   % Find the elements that are coupled at least as much as this increment (with a small correction)
   coupinc = coup >= incr(i)-2*(id/10);
   incr(i)
   % Determine how many clusters comprise this increment
   clus = meshclusters(p, coupinc);
   % Loop through clusters to define individual events
   for j = 1:length(clus)
      % Multiply binary array by element area and sum to give total area
      tarea(j, i) = sum(area(clus{j}));
   end
end
% Trim unused rows
tarea(sum(isnan(tarea), 2) == length(incr), :) = [];

% Multiply total area by Wells and Coppersmith (1994) scaling parameters to give M_W

% Check to see if Wells and Coppersmith scaling parameters were specified as inputs
if ~exist('wcscale', 'var')
   % If not, set to be equal to the "all" coefficients
   a = 4.07;
   b = 0.98;
   aunc = 0.06;
   bunc = 0.03;
else
   % If so,
   
   % If it's numeric, then values are specified
   if isnumeric(wcscale)
      [a, b, aunc, bunc] = deal(wcscale(1), wcscale(2), wcscale(3), wcscale(4));
   elseif ischar(wcscale)
      switch lower(wcscale)
         case 'all'
            [a, b, aunc, bunc] = deal(4.07, 0.98, 0.06, 0.03);
         case 'n'
            [a, b, aunc, bunc] = deal(3.93, 1.02, 0.23, 0.10);
         case 'r'
            [a, b, aunc, bunc] = deal(4.33, 0.90, 0.12, 0.05);
         case 'ss'
            [a, b, aunc, bunc] = deal(3.98, 1.02, 0.07, 0.03);
      end
   end
end
% Convert area to moment magnitude
mmag = a + b*log10(tarea);
mmag(mmag < 7.0) = NaN;
merr = aunc + bunc*log10(tarea);

% Make a plot
figure
% Define min. and max. range for each increment value
minv = min(mmag-merr, [], 1);
maxv = max(mmag+merr, [], 1);
range = patch([incr(:); flipud(incr(:))], [minv(:); flipud(maxv(:))], 0.75*[1 1 1]);
set(range, 'edgecolor', 'none');
hold on
if ~exist('subset', 'var')
   subset = 1:length(incr);
end
plot(incr(:, subset), mmag(:, subset), '.w', 'markersize', 25, 'markeredgecolor', 'k');
axis([minmax(incr) 7.0 9.5]);
set(gca, 'xticklabel', num2str((0.1:0.1:1)', '%.1f'))
set(gca, 'yticklabel', num2str((7:0.5:9.5)', '%.1f'))
prepfigprint
xlabel('Coupling fraction', 'fontsize', 12)
ylabel('M_W', 'fontsize', 12)


