function matinv = openmeeg_inv(mat)

% FT_HEADMODEL_OPENMEEG creates a volume conduction model of the
% head using the boundary element method (BEM). This function takes
% as input the triangulated surfaces that describe the boundaries and
% returns as output a volume conduction model which can be used to
% compute leadfields.
%
% This function implements 
%   Gramfort et al. OpenMEEG: opensource software for quasistatic
%   bioelectromagnetics. Biomedical engineering online (2010) vol. 9 (1) pp. 45
%   http://www.biomedical-engineering-online.com/content/9/1/45
%   doi:10.1186/1475-925X-9-45
% and
%   Kybic et al. Generalized head models for MEG/EEG: boundary element method
%   beyond nested volumes. Phys. Med. Biol. (2006) vol. 51 pp. 1333-1346
%   doi:10.1088/0031-9155/51/5/021
% 
% The implementation in this function is derived from the the OpenMEEG project
%  and uses external command-line executables. See http://gforge.inria.fr/projects/openmeeg
% and http://gforge.inria.fr/frs/?group_id=435.
%
% Use as
%   vol = ft_headmodel_openmeeg(geom, ...)
%
% Optional input arguments should be specified in key-value pairs and can
% include
%   hdmfile          = string, filename with BEM headmodel
%   conductivity     = vector, conductivity of each compartment
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

%$Id: ft_headmodel_openmeeg.m 8305 2013-07-02 09:59:57Z roboos $

ft_hastoolbox('openmeeg', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this uses an implementation that was contributed by INRIA Odyssee Team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% show the license once
% openmeeg_license

% check that the binaries are ok
om_checkombin;

% store the current path and change folder to the temporary one
tmpfolder = cd;

try
  cd(tempdir)
 
  % these will hold the shell script and the inverted system matrix
  [tmp,tname] = fileparts(tempname);
  if ~ispc
    exefile = [tname '.sh'];
  else
    exefile = [tname '.bat'];
  end
  
  [tmp,tname] = fileparts(tempname);
  hmfile    = [tname '.bin'];
  [tmp,tname] = fileparts(tempname);
  hminvfile = [tname '.bin'];
  
  % write conductivity and geometry files
  om_save_sym(mat,hmfile);
  
  % Exe file
  efid = fopen(exefile, 'w');
  omp_num_threads = feature('numCores');
  if ~ispc
    fprintf(efid,'#!/usr/bin/env bash\n');
    fprintf(efid,['export OMP_NUM_THREADS=',num2str(omp_num_threads),'\n']);
    fprintf(efid,['om_minverser ./' hmfile ' ./' hminvfile ' 2>&1 > /dev/null\n']);
  else
    fprintf(efid,['om_minverser ./' hmfile ' ./' hminvfile '\n']);
  end
  
  fclose(efid);
  
  if ~ispc
    dos(sprintf('chmod +x %s', exefile));
  end
catch
  cd(tmpfolder)
  rethrow(lasterror)
end

try
  % execute OpenMEEG and read the resulting file
  if ispc
    dos([exefile]);
  else
    version = om_getgccversion;
    if version>3
      dos(['./' exefile]);
    else
      error('non suitable GCC compiler version (must be superior to gcc3)');
    end
  end
  matinv = om_load_sym(hminvfile,'binary');
  cleaner(hmfile,hminvfile,exefile)
  cd(tmpfolder)
catch
  warning('an error ocurred while running OpenMEEG');
  disp(lasterr);
  cleaner(vol,bndfile,condfile,geomfile,hmfile,hminvfile,exefile)
  cd(tmpfolder)
end

% remember the type of volume conduction model
vol.type = 'openmeeg';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cleaner(hmfile,hminvfile,exefile)
% delete the temporary files
delete(hmfile);
delete(hminvfile);
delete(exefile);

function ok = checknormals(bnd)
% FIXME: this method is rigorous only for star shaped surfaces
ok = 0;
pnt = bnd.pnt;
tri = bnd.tri;
% translate to the center
org = mean(pnt,1);
pnt(:,1) = pnt(:,1) - org(1);
pnt(:,2) = pnt(:,2) - org(2);
pnt(:,3) = pnt(:,3) - org(3);

w = sum(solid_angle(pnt, tri));

if w<0 && (abs(w)-4*pi)<1000*eps
  ok = 0;
  warning('your normals are outwards oriented\n')
elseif w>0 && (abs(w)-4*pi)<1000*eps
  ok = 1;
%   warning('your normals are inwards oriented')
else
  error('your surface probably is irregular\n')
  ok = 0;
end
