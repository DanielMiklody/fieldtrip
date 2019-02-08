function om_write_geom(geomfile,bndfile,names)
%   OM_WRITE_COND
%       [] = OM_WRITE_GEOM(GEOMFILE,BNDFILE)
%       [] = OM_WRITE_GEOM(GEOMFILE,BNDFILE,NAMES)
%
%   Write geometry file for OpenMEEG
%
%   Authors: Alexandre Gramfort alexandre.gramfort@inria.fr
%            Paul Czienskowski

% Copyright (C) 2010-2017, OpenMEEG developers

numbnd = length(bndfile);

if nargin < 3
    names = {};
    for k=1:numbnd
        names{k} = ['domain', num2str(k)];
    end
    
    gfid = fopen(geomfile, 'w');
    
    if gfid == -1
        error(['Failed to open file ''',geomfile,'''.']);
    end
    
    fprintf(gfid,'# Domain Description 1.0\n');
    fprintf(gfid,'                        \n');
    fprintf(gfid,'Interfaces %d Mesh      \n', numbnd);
    fprintf(gfid,'                        \n');
    
    for i=1:numbnd
        fprintf(gfid,'%s                  \n', bndfile{i});
    end
    
    fprintf(gfid,'                        \n');
    fprintf(gfid,'Domains %d              \n', numbnd+1);
    fprintf(gfid,'                        \n');
    
    fprintf(gfid,'Domain air %d           \n', 1);
    
    
    for i=1:numbnd
        if i < numbnd
            fprintf(gfid,'Domain %s %d %d\n', names{i}, i+1, -i);
        else
            fprintf(gfid,'Domain %s %d   \n', names{i}, -i);
        end
    end
    fclose(gfid);
else
    if numbnd ~= length(names)
        error('Number of boundary files is not equal to the number of domain names.');
    end
    gfid = fopen(geomfile, 'w');
    
    if gfid == -1
        error(['Failed to open file ''',geomfile,'''.']);
    end
    %     Inames=strcat('I',names);
    %     Inamesp=strcat('+I',names);
    %     Inamesn=strcat('-I',names);
    for ii=1:numbnd
            Numbernames{ii}=num2str(ii);
        switch names{ii}
            case {'skin','scalp'}
                indices.skin=ii;
                %                 names{ii}='skin';
            case {'skull','bone'}
                %                 names{ii}='skull';
                indices.skull=ii;
            case {'eyel','left eye'}
                %                 names{ii}='eyel';
                indices.eyel=ii;
            case {'eyer','right eye'}
                %                 names{ii}='eyer';
                indices.eyer=ii;
            case {'csf','CSF'}
                %                 names{ii}='csf';
                indices.csf=ii;
            case {'brain','gray matter'}
                %                 names{ii}='brain';
                indices.gray=ii;
            case {'white matter'}
                %                 names{ii}='brain';
                indices.white=ii;
        end
    end
    Inamesp=strcat('+',Numbernames);
    Inamesn=strcat('-',Numbernames);
    if ~isfield(indices,'csf')
        names{end+1}='';
        Inamesp{end+1}=['+I' names{indices.gray}];
        Inamesn{end+1}=['-I' names{indices.gray}];
        indices.csf=numel(names);
    end
    if ~isfield(indices,'white')
        names{end+1}='';
        Inamesp{end+1}='';
        Inamesn{end+1}='';
        indices.white=numel(names);
    end
    if ~isfield(indices,'eyer')
        names{end+1}='';
        Inamesp{end+1}='';
        Inamesn{end+1}='';
        indices.eyer=numel(names);
    end
    if ~isfield(indices,'eyel')
        names{end+1}='';
        Inamesp{end+1}='';
        Inamesn{end+1}='';
        indices.eyel=numel(names);
    end
    
    fprintf(gfid,'# Domain Description 1.0\n');
    fprintf(gfid,'                        \n');
    fprintf(gfid,'Interfaces %d Mesh      \n', numbnd);
    fprintf(gfid,'                        \n');
    
    for i=1:numbnd
        %         fprintf(gfid,'Interface: "%s"                  \n',Inames{i}, bndfile{i});
        fprintf(gfid,'%s                  \n', bndfile{i});
    end
    
    fprintf(gfid,'                        \n');
    fprintf(gfid,'Domains %d              \n', numbnd+1);
    fprintf(gfid,'                        \n');
    
    fprintf(gfid,'Domain air %s           \n', Inamesp{indices.skin});
    for ii=1:numbnd
        switch names{ii}
            case {'skin','scalp'}
                fprintf(gfid,'Domain %s %s %s %s %s\n', names{ii}, Inamesp{indices.skull},  Inamesn{indices.skin},  Inamesp{indices.eyel},  Inamesp{indices.eyer});
            case {'skull','bone'}
                fprintf(gfid,'Domain %s %s %s\n', names{ii}, Inamesp{indices.csf}, Inamesn{indices.skull});
            case {'eyel','left eye'}
                fprintf(gfid,'Domain %s %s   \n', names{ii}, Inamesn{indices.eyel});
            case {'eyer','right eye'}
                fprintf(gfid,'Domain %s %s   \n', names{ii}, Inamesn{indices.eyer});
            case {'csf','CSF'}
                fprintf(gfid,'Domain %s %s %s\n', names{ii}, Inamesp{indices.gray}, Inamesn{indices.csf});
            case {'brain','gray matter'}
                fprintf(gfid,'Domain %s %s %s \n', names{ii}, Inamesn{indices.gray}, Inamesp{indices.white});
            case {'white matter'}
                fprintf(gfid,'Domain %s %s   \n', names{ii}, Inamesn{indices.white});
        end
    end
    fclose(gfid);
end
end %  function