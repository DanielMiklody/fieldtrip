function [gain] = openmeeg_gain(elec, vol)

% OPENMEEG_EIT computes the OpenMEEG EIT matrix
%              i.e. Right hand side in the potential equation
%
% Use as
%   [eit] = openmeeg_eit(po, vol, flag)
%
% flag = 1 non adaptive algorithm: does not try to approximate the
% potential in the neighborhodd of the leads, by locally refining the BEM surface

% Copyright (C) 2009, Alexandre Gramfort
% INRIA Odyssee Project Team

% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailed information


% store the current path and change folder to the temporary one
tmpfolder = cd;
om_checkombin;

bndom = vol.bnd;

try
    cd(tempdir)

    % write the triangulations to file
    bndfile = {};
    for i=1:length(vol.bnd)
        [~,tname] = fileparts(tempname);
        bndfile{i} = [tname '.tri'];
        ok = checknormals(bndom(i));
        if ~ok
          bndom(i).tri = fliplr(bndom(i).tri);
        end
        om_save_tri(bndfile{i}, bndom(i).pnt, bndom(i).tri);
    end
    
    % these will hold the shell script and the inverted system matrix
    [~,tname] = fileparts(tempname);
    if ~ispc
      exefile = [tname '.sh'];
    else
      exefile = [tname '.bat'];
    end

    [~,tname] = fileparts(tempname);
    condfile = [tname '.cond'];
    [~,tname] = fileparts(tempname);
    geomfile = [tname '.geom'];
    [~,tname] = fileparts(tempname);
    electrodefile = [tname '.elec'];
    [~,tname] = fileparts(tempname);
    eitsourcefile = [tname '.bin'];
    [~,tname] = fileparts(tempname);
    gainfile = [tname '.bin'];
    [~,tname] = fileparts(tempname);
    hmvinvfile = [tname '.bin'];
    [~,tname] = fileparts(tempname);
    h2mfile = [tname '.bin'];

    % write conductivity and geometry files
    om_write_geom(geomfile,bndfile);
    om_write_cond(condfile,vol.cond);
    om_save_full(eye(size(vol.mat)),eitsourcefile);
    om_save_sym(eye(size(vol.mat)),hmvinvfile);

    % handle electrode file
    elec=ft_convert_units(elec,vol.unit);
    om_save_full(elec.elecpos,electrodefile,'ascii');

    % Exe file
    efid = fopen(exefile, 'w');
    omp_num_threads = feature('numCores');

    str2 = ' -H2EM';
    str3 = ' -EEG';

    if ~ispc
      fprintf(efid,'#!/usr/bin/env bash\n');
      fprintf(efid,['export OMP_NUM_THREADS=',num2str(omp_num_threads),'\n']);
      fprintf(efid,['om_assemble' str2 ' ./',geomfile,' ./',condfile,' ./',electrodefile,' ./',h2mfile,' 2>&1 > /dev/null\n']);
      fprintf(efid,['om_gain' str3 ' ./',hmvinvfile,' ./',eitsourcefile,' ./',gainfile,' 2>&1 > /dev/null\n']);
    else
      fprintf(efid,['om_assemble' str2 ' ./',geomfile,' ./',condfile,' ./',electrodefile,' ./',h2mfile,'\n']);
      fprintf(efid,['om_gain' str3 ' ./',hmvinvfile,' ./',eitsourcefile,' ./',h2mfile,' ./',gainfile,' \n']);
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
    disp(['Assembling OpenMEEG Sensor Interpolation matrix']);
    stopwatch = tic;
    if ispc
        dos([exefile]);
    else
        dos(['./' exefile]);
    end
    gain = om_load_full(gainfile,'binary');
    toc(stopwatch);
    cleaner(vol,bndfile,condfile,geomfile,exefile,electrodefile,eitsourcefile,hmvinvfile)
    cd(tmpfolder)
catch
    warning('an error ocurred while running OpenMEEG');
    disp(lasterr);
    cleaner(vol,bndfile,condfile,geomfile,exefile,electrodefile,eitsourcefile,hmvinvfile)
    cd(tmpfolder)
end

function cleaner(vol,bndfile,condfile,geomfile,exefile,electrodefile,eitsourcefile,hmvinvfile)
% delete the temporary files
for i=1:length(vol.bnd)
    if exist(bndfile{i},'file'),delete(bndfile{i}),end
end
if exist(condfile,'file'),delete(condfile);end
if exist(geomfile,'file'),delete(geomfile);end
if exist(exefile,'file'),delete(exefile);end
if exist(electrodefile,'file'),delete(electrodefile);end
if exist(eitsourcefile,'file'),delete(eitsourcefile);end
if exist(hmvinvfile,'file'),delete(hmvinvfile);end

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
%   warning('your normals are outwards oriented\n')
elseif w>0 && (abs(w)-4*pi)<1000*eps
  ok = 1;
%   warning('your normals are inwards oriented')
else
  error('your surface probably is irregular\n')
  ok = 0;
end