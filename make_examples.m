function make_examples(dirs,filename,str)
%MAKE_EXAMPLES  Publish the files in the examples directory to website.
% MAKE_EXAMPLES(DIR,FILENAME,'NEW') will only publish the file FILENAME.M in
% directory $Chebfunroot/examples/DIR/ where FILENAME and DIR must be
% strings. The additional string 'NEW' creates the additional index files,
% tags, and listings that are required for new Examples.

% Admittedly, this file is a mess and should be removed from releases (or
% at least replaced with a more user friendly code that simply makes html
% versions of the Examples).

% The flags below can only be adjusted manually.
html = false;  % Publish to html? (This should be true when released).
pdf = false;   % By default this will be off.
shtml = true; % This should only be used by admin for creating the 
              % shtml files for the Chebfun website.
clean = false;
listing = false;
exampleshtml = false;
release = false;

nargs = nargin;

webdir = '/common/htdocs/www/maintainers/hale/chebfun/';
% if ~isunix || ~exist(webdir,'dir')
%     error(['This file is designed only to run on a Linux box ' ...
%         'which has access to the Oxford webserver.']);
% end

currUser = getenv('USER');
adminUser = 'hale';
if strcmp(currUser,adminUser), admin = true; else admin = false; end

precsign = '%';

if nargin > 0 && ischar(dirs) 
    if strncmp(dirs,'listing',4); listing = true; end
    if strcmp(dirs,'clean'); clean = true; end
    if strcmp(dirs,'tags'); make_tags(); return, end
    
    if strcmp(dirs,'permissions'); 
        cd(webdir)
        eval('!chgrp -Rf chebfun *')
        eval('!chmod -Rf 775 *')
        eval('!chown -Rf hale *')
        return
    end
   
    if strcmp(dirs,'html'); 
        exampleshtml = true; 
        shtml = false;
        html = true;
        pdf = false;
        listing = false;
    end
    if strcmp(dirs,'release'); 
        release = true;
        exampleshtml = false; 
        shtml = false;
        html = true;
        pdf = false;
        listing = false;
    end
    if nargin > 1 && strcmp(filename,'make'); 
        exampleshtml = false; 
        shtml = true;
        html = true;
        pdf = true;
        listing = false;
        nargs = 1;
    end
    if nargin == 3 && strcmpi(str,'new')
        make_examples(dirs,filename);
        make_examples(dirs);
        make_examples('list');
        make_tags();
%         make_examples('permissions');
        return
    end        
end

examplesdir = pwd;
% Define formatting.
% HTML formatting
opts = [];
opts.stylesheet = fullfile(examplesdir,'templates','custom_mxdom2simplehtml.xsl');
if ~shtml || ~exist(opts.stylesheet,'file'), 
    opts.stylesheet = [];   % Resort to default stylesheet.
end       
opts.catchError = false;
% opts.evalCode = false;

% PDF formatting
optsPDF = [];
optsPDF.stylesheet = fullfile(examplesdir,'templates','custom_mxdom2latex.xsl'); 
% optsPDF.format = 'pdf'; 
optsPDF.format = 'latex'; 
if ~shtml || ~exist(optsPDF.stylesheet,'file') || strcmp(optsPDF.format,'pdf'), 
    optsPDF.stylesheet = [];   % Resort to default stylesheet.
end 
optsPDF.outputDir = 'pdf'; 
optsPDF.catchError = false;

% These examples are special because they produce output from the
% anon/display field, which, because java is running will attempt to pipe
% hyperlinks to the command window. This results in a massive mess, which
% we clean up by filtering the html and tex files.
javalist = {'ChebfunAD'};

if nargs == 0 || listing || clean || exampleshtml || release
    % Find all the directories (except some exclusions).
    dirlist = struct2cell(dir(fullfile(pwd)));
    dirs = {}; k = 1;
    for j = 3:size(dirlist,2)
        if dirlist{4,j} && ...
                ~strncmp(dirlist{1,j},'.',1) && ...
                ~strcmp(dirlist{1,j},'templates') && ...
                ~strcmp(dirlist{1,j},'old_examples') && ...
                ~strcmp(dirlist{1,j},'old') && ...
                ~strcmp(dirlist{1,j},'temp')
            dirs{k} = dirlist{1,j};
            k = k + 1;
        end
    end
elseif nargs == 1 && ischar(dirs)
    % Directory has been passed
    dirs = {dirs};
elseif nargs == 2
    html = true;
    pdf = true;
    % Compile a single file (given)
    if iscell(dirs), dirs = dirs{:}; end
    if strcmp(filename(end-1:end),'.m'), 
        filename = filename(1:end-2); 
    end
    cd(dirs);
    %%% HTML %%%
    if html
        % Clean up
        if exist('html','dir')
            delete(['html/',filename,'*.png'])
        end
        
        % Publish (to dirname/html/dirname.html)
        mypublish([filename,'.m'],opts); 
        
        % Make the filename clickable
        cd html
        curfile = [dirs,'/',filename,'.m']; 
        filetext = fileread([filename,'.html']);
        if shtml
            newtext = sprintf('<a href="/chebfun/examples/%s">%s</a>', ...
                curfile,curfile);
        else
            newtext = sprintf('<a href="%s">%s</a>', ...
                fullfile(examplesdir,dirs,[filename,'.m']),curfile);
        end
        filetext = strrep(filetext,curfile,newtext);
        
        % Deal with hashtags
        filetext = strrep(filetext,[newtext ') [Tags:'],[newtext ')<br/>[Tags:']);
        idx1 = strfind(lower(filetext),'[tags:');
        if ~isempty(idx1)
            idx2 = strfind(filetext(idx1:end),'#');
            idx3 = strfind(filetext(idx1:end),']');
            tagline = filetext(idx1(1)-1+(idx2(1):idx3(1)-1));
            idx = strfind(tagline,'#');
        else
            idx = [];
        end
        if ~isempty(idx)
            idx2 = strfind(tagline,',');
            tags = {};
            for k = 1:numel(idx2)
                tags{k} = tagline(idx(k)+1:idx2(k)-1);
            end
            tags{numel(idx)} = tagline(idx(end)+1:end);
            newtagline = [];
            for k=1:numel(tags)
                url = '/chebfun/examples/tags.php?query=';
                newtagline = [newtagline '<a href="' url tags{k} '">#' tags{k} '</a>, '];
            end
            newtagline = newtagline(1:end-2);
            filetext = strrep(filetext,tagline,newtagline);
        end   
        
%         % Copyright notice
%         if shtml
%             filetext = strrep(filetext,'<p>Licensed under a Creative Commons 3.0 Attribution license','');
%             filetext = strrep(filetext,'<a href="http://creativecommons.org/licenses/by/3.0/">http://creativecommons.org/licenses/by/3.0/</a>','');
%             filetext = strrep(filetext,'by the author above.</p>','');
%         end

        if any(strcmp(filename,javalist))
%             filetext = strrep(filetext,'&lt;','<');
%             filetext = strrep(filetext,'&gt;','>');
                starts = strfind(filetext,'%&lt;a href="matlab: edit');
                for k = numel(starts):-1:1
                    endsk = strfind(filetext(starts(k):starts(k)+100),'&lt;/a');
                    strk = filetext(starts(k)+(0:endsk(1)+2));
                    idx1 = strfind(strk,'&gt;'); idx2 = strfind(strk,'&lt;');
                    newstr = strk(idx1(1)+4:idx2(2)-1);
                    filetext(starts(k)+(1:numel(newstr))) = newstr;
                    filetext(starts(k)+((numel(newstr)+1):endsk(1)+8)) = [];
%                         filetext(starts(k)+(0:endsk(1)+2)) = newstr;
                end
        end
        
        if strcmp(filename,'Writing3D')
            png = 'Writing3D_04.png';
            gif = 'Writing3D_04.gif';
            filetext = strrep(filetext,png,gif);
        end
        
        if strcmp(filename,'BouncingBall')
            png = 'BouncingBall_01.png';
            gif = 'BouncingBall_01.gif';
            try
                eval(['!cp BouncingBall.gif ' webdir 'examples/ode/html/BouncingBall_01.gif'])
            end
            filetext = strrep(filetext,png,gif);
        end
        
        if strcmp(filename,'InteractiveInterp')
            png = 'InteractiveInterp_01.png';
            gif = 'InteractiveInterp_01.gif';
            try
                eval(['!cp ../InteractiveInterp_01.gif ' webdir 'examples/approx/html/'])
            end
            filetext = strrep(filetext,png,gif);
        end
        
        if strcmp(filename,'Scribble2')
            png = 'Scribble2_03.png';
            gif = 'Scribble2_03.gif';
            try
                eval(['!cp ../Scribble2.gif ' webdir 'examples/fun/html/Scribble2_03.gif'])
            end
            filetext = strrep(filetext,png,gif);
        end
        
        % Try to insert hyperlinks for references
        try
            refloc = strfind(lower(filetext),'reference');
            sourceloc = strfind(filetext,'SOURCE BEGIN');
            if isempty(sourceloc)
                reftext = filetext(refloc:end);
                sourcetext = [];
            else
                refloc(refloc>sourceloc) = [];
                refloc = refloc(end);
                reftext = filetext(refloc:sourceloc);
                sourcetext = filetext(sourceloc+1:end);
            end
            linklocs1 = strfind(reftext,'<a href'); 
            linklocs2 = strfind(reftext,'>');
            k = 1;
            while 1
                refk = ['[' int2str(k) ']'];
                idxk = strfind(reftext,refk);
                if isempty(idxk), break, end
                idx(k) = idxk;
                k = k+1;
            end
            bodytext = filetext(1:refloc-1);
            idx(k) = inf;
            for k = 1:numel(idx)-1
                startk = linklocs1(linklocs1>idx(k) & linklocs1<idx(k+1));
                endk = linklocs2(linklocs2>idx(k) & linklocs2<idx(k+1));
                if isempty(startk) || isempty(endk), continue, end
                strk = reftext(startk:endk);
                refk = ['[' int2str(k) ']'];
                bodytext = strrep(bodytext,refk,[strk,refk,'</a>']);
            end
            filetext = [bodytext reftext sourcetext];
        end

        fidhtml = fopen([filename,'.html'],'w+');
        fprintf(fidhtml,'%s',filetext);
        fclose(fidhtml);
        cd ..
    end
    %%% PDF %%%
    if pdf
        % Publish (to dirname/pdf/dirname.pdf)
        mypublish([filename,'.m'],optsPDF);        
        if strcmp(optsPDF.format,'latex') && isunix
            try
                cd pdf
                
                % Tidy up special characters
                filetext = fileread([filename,'.tex']);
                filetext = strrep(filetext,'ü','\"{u}');
                filetext = strrep(filetext,'ø','{\o}');
                filetext = strrep(filetext,'ó','\''{o}');
                filetext = strrep(filetext,'ö','\"{o}');
                filetext = strrep(filetext,'ő','\H{o}');
                filetext = strrep(filetext,'Ő','\H{O}');  
                filetext = strrep(filetext,'é','\''{e}');  
                
                % Fix a MATLAB bug!
                filetext = strrep(filetext,'$\$','$$');  

                if any(strcmp(filename,javalist))
                    starts = strfind(filetext,'%<a href="matlab: edit');
                    for k = numel(starts):-1:1
                        endsk = strfind(filetext(starts(k):starts(k)+100),'</a>');
                        strk = filetext(starts(k)+(0:endsk(1)+2));
                        idx1 = strfind(strk,'>'); idx2 = strfind(strk,'<');
                        newstr = strk(idx1(1)+1:idx2(2)-1);
                        filetext(starts(k)+(1:numel(newstr))) = newstr;
                        filetext(starts(k)+((numel(newstr)+1):endsk(1)+2)) = [];;
%                         filetext(starts(k)+(0:endsk(1)+2)) = newstr;
                    end
                end
                
                fidpdf = fopen([filename,'.tex'],'w+');
                fprintf(fidpdf,'%s',filetext);
                fclose(fidpdf);
                eval(['!latex ', filename, '> foo.tmp'])
                eval(['!rm foo.tmp'])
                eval(['!dvipdfm -q ', filename])                
%                 ! rm *.aux *.log *.tex  *.dvi
                cd ../
            catch
                warning('CHEBFUN:examples:PDFfail','PDF PUBLISH FAILED.');
            end            
        end
    end
    cd ..
    
    % Upload to web server.
    if shtml
        curdir = pwd;
        cd(fullfile(webdir,'examples'))
        if ~exist(dirs,'dir'), mkdir(dirs), end
        cd(dirs)
        if ~exist('html','dir'), mkdir('html'), end
        if ~exist('pdf','dir'), mkdir('pdf'), end
        
        fprintf('Uploading html...')
        eval(['!cp ', fullfile(curdir,dirs,'html',[filename,'.html']),' ','html'])
        try
            eval(['!cp ', fullfile(curdir,dirs,'html',[filename,'*.png']),' ','html'])
        end
        if exist(fullfile(curdir,dirs,'html',[filename,'.shtml']),'file')
            eval(['!cp ', fullfile(curdir,dirs,'html',[filename,'.shtml']),' ','html'])
        end
        fprintf('Complete. ')
        
        fprintf('Uploading pdf...')
        eval(['!cp ',fullfile(curdir,dirs,'pdf',[filename,'.pdf']),' ',fullfile('pdf',[filename,'.pdf'])])
        fprintf('Complete. ')
        
        fprintf('Uploading m files...')
        eval(['!cp ',fullfile(curdir,dirs,[filename,'.m']),' ',[filename,'.m']]);
        fprintf('Complete. ')
        
        fprintf('Setting file permissions...')
            eval(['!chgrp -f chebfun ' filename '.m'])
            eval(['!chmod -f 775 ' filename '.m'])
            eval(['!chown -f hale ' filename '.m'])
        cd html
            eval(['!chgrp -f chebfun ' filename '*'])
            eval(['!chmod -f 775 ' filename '*'])
            eval(['!chown -f hale ' filename '*'])
        cd ../pdf/
            eval(['!chgrp -f chebfun ' filename '.pdf'])
            eval(['!chmod -f 775 ' filename '.pdf'])
            eval(['!chown -f hale ' filename '.pdf'])
        fprintf('Complete.\n')
        
        cd(curdir)
    end
    return
end

% Clean up
if clean
    fprintf('Cleaning. Please wait ...\n')
    for j = 1:numel(dirs)
        if ~exist(dirs{j},'dir'), continue, end
        cd(dirs{j})
        cd
        delete *.html *.shtml
        if exist('html','dir'), rmdir('html','s'), end
        if exist('pdf','dir'), rmdir('pdf','s'), end
        cd ..
    end
    fprintf('Done.\n')
    return
end

% Make index
if listing
    fprintf('Compiling index. Please wait ...\n')
    
    mfile = {}; filedir = {};
    % Find *all* the files.
    for j = 1:numel(dirs)
        % Move to the directory.
        if strcmp(dirs{j},'temp'), continue, end
        cd(dirs{j})
        % Find all the m-files.
        dirlist = dir(fullfile(pwd,'*.m'));
        mfile = [mfile dirlist.name];
        filedir = [filedir cellstr(repmat(dirs{j},numel(dirlist),1))'];
        cd ..
    end

    % Get the ordering (we need to ignore A and THE)
    desc = cell(numel(mfile),1); 
    for k = 1:numel(mfile)
        filename = mfile{k}(1:end-2);
        cd(filedir{k})
        % Grab the file description.
        fidk = fopen([filename,'.m']);
        txt = fgetl(fidk);
        origtxt = txt;
        txt = upper(txt);
        fclose(fidk);
        if txt < 1, continue, end % This mfile will be ignored
        if numel(txt) >1 && strcmp(txt(1:2),'%%')
            txt = txt(4:end);
        else
            txt = '     ';  % This mfile will be ignored
        end
        if strcmpi(txt(1:2),'A ')
            txt = txt(3:end);
        elseif strcmpi(txt(1:3),'AN ')
            txt = txt(4:end);             
        elseif strcmpi(txt(1:4),'THE ')
            txt = txt(5:end);
        end
        desc{k} = txt;
        origtxt = origtxt(4:end);
        if numel(origtxt) > 50
            idx = strfind(origtxt,':');
            if ~isempty(idx), origtxt = origtxt(1:idx(1)-1); end
        end
        origdesc{k} = origtxt;
        cd ..
    end
    [desc indx] = sort(desc);
    origdesc = origdesc(indx);
    mfile = mfile(indx);
    filedir = filedir(indx);
    
%     % Print data to file (list version)
%     fid = fopen('listing.html','w');
%     fprintf(fid,'<ul class="atap" style="padding-left:15px;">\n');
%     for k = 1:numel(mfile)
%          if isempty(desc{k}), continue, end
%          if strcmp(desc{k}(1),' '), continue, end
%          origdesc{k} = capitalize(origdesc{k});
%          mfile{k} = mfile{k}(1:end-2);
%          newtext = sprintf(['  <li>%s  <span style="float:right"><a href="%s/" ',...
%              'style="width:70px; display: inline-block;">%s</a>', ...
%              '(<a href="%s/html/%s.shtml">html</a>, <a href="%s/pdf/%s.pdf">PDF</a>, ',...
%              '<a href="%s/%s.m">M-file</a>)</span></li>\n\n'], ...
%                     origdesc{k},filedir{k},filedir{k},filedir{k},mfile{k},filedir{k},mfile{k},filedir{k},mfile{k});
%          fprintf(fid,newtext);
%     end
%     fclose(fid);
    
    % Print data to file (table version)
    fid = fopen('listing.html','w');
    if ~shtml
        fprintf(fid,'Here is the complete list of Chebfun Examples and the sections they belong to.<br/><br/>\n');
    end
    fprintf(fid,'<table style="padding-left:15px; cellpadding:2px; width:100%s;">\n',precsign);
    ms = 0;
    for k = 1:numel(mfile)
         if isempty(desc{k}), continue, end
         if strcmp(desc{k}(1),' '), continue, end
%          origdesc{k} = capitalize(origdesc{k});
         mfile{k} = mfile{k}(1:end-2);
         ms = max(ms,length(origdesc{k}));
         if shtml
%              newtext = sprintf(['  <tr>\n   <td style="text-transform: uppercase;">%s</td>\n', ...
             newtext = sprintf(['  <tr>\n   <td>%s</td>\n', ...
                 '   <td style="float:right"><a href="%s/" style="width:70px; display: inline-block;">%s</a></td>\n', ...
                 '   <td>(<a href="%s/html/%s.shtml">html</a>, <a href="%s/pdf/%s.pdf">PDF</a>, ',...
                 '<a href="%s/%s.m">M-file</a>)</td>\n  </tr>\n\n'], ...
                        origdesc{k},filedir{k},filedir{k},filedir{k},mfile{k},filedir{k},mfile{k},filedir{k},mfile{k});
         else
%              newtext = sprintf(['  <tr>\n   <td style="text-transform: uppercase;"><a href="%s/%s.m">%s</a></td>\n', ...
             newtext = sprintf(['  <tr>\n   <td><a href="%s/%s.m">%s</a></td>\n', ...
                 '   <td style="float:right">(<a href="%s/" style="width:70px; display: inline-block;">%s</a>)</td>\n  </tr>\n\n'], ...
                        filedir{k},mfile{k},origdesc{k},filedir{k},filedir{k});
         end
         fprintf(fid,newtext);
    end
    fprintf(fid,'<table>\n');
    fclose(fid);
    
%     % Print data to file (.txt version)
%     fid = fopen('LIST.txt','w');
%     fprintf(fid,'Below is a list of all the available Chebfun Examples in this directory\nMore can be found on the web at http://www.maths.ox.ac.uk/chebfun/examples/\n\n');
%     for k = 1:numel(mfile)
%          if isempty(desc{k}), continue, end
%          if strcmp(desc{k}(1),' '), continue, end
%          origdesc{k} = capitalize(origdesc{k});
%          ws = repmat(' ',1,ms-length(origdesc{k})+4);
% %          mfile{k} = mfile{k}(1:end-2);
%          newtext = sprintf('%s%s%s/%s.m\n',origdesc{k},ws,filedir{k},mfile{k});
%          fprintf(fid,newtext);
%     end
%     fclose(fid);

    if shtml 
        curdir = pwd;
        cd(webdir)
        fprintf([' Uploading. Please wait ... '])
        eval(['!cp ',fullfile(curdir,'listing.html'),' ','examples']);
        eval(['!chgrp -f chebfun examples/listing.html']);
        eval(['!chmod -f 775 examples/listing.html']);
        eval(['!chown -f hale examples/listing.html']);
        cd(curdir)
        fprintf('Done.\n')
    end
    
    fprintf('Done.\n')
    return
end

% % Make examples.html
fid0 = fopen('examples.html','w+');
% Open template
fid_et1 = fopen('templates/examples_template1.txt','r');
% Read data.
tmp = fread(fid_et1,inf,'*char');
fclose(fid_et1);
% Write
fprintf(fid0,' %s',tmp);

% Sort the directories to match contents.txt
if numel(dirs) > 1 %|| iscell(dirs) && numel(dirs{1})
    fidc = fopen('contents.txt','r+');
    titletxt = fgetl(fidc);
    titles = [];
    while titletxt > 0
        titles = [titles ; titletxt(1:3)];
        titletxt = fgetl(fidc);
    end
    titles = cellstr(titles);
    titles(strcmp(titles,'tem')) = [];
    if numel(dirs) == numel(titles)
        [ignored idx1] = sort(titles);
        [ignored idx2] = sort(idx1);
        dirs = dirs(idx2);
    end        
end

if exampleshtml
    for j = 1:numel(dirs) 
        % Find the title of this directory
        fidc = fopen('contents.txt','r+');
        titletxt = fgetl(fidc);
        while ~strncmp(dirs{j},titletxt,3)
            titletxt = fgetl(fidc);
            if titletxt < 0, 
                error('CHEBFUN:examples:dirname', ...
                    ['Unknown directory name "',dirs{j},'. Update contents.txt.']);
            end        
        end
        titletxt = titletxt(length(dirs{j})+2:end);
        % Add entry to examples/examples.html
        if ~strcmp(dirs{j},'temp')
            fprintf(fid0,['<li><a href="',dirs{j},'/',dirs{j},'.html" style="text-transform: uppercase;">',titletxt,'</a>\n</li>\n\n']);
        end
    end
    % Open template
    fid_et2 = fopen('templates/examples_template2.txt','r');
    % Read data.
    tmp = fread(fid_et2,inf,'*char');
    fclose(fid_et2);
    % Write
    fprintf(fid0,' %s',tmp);
    fclose(fid0);
    return
end

% Loop over the directories.
for j = 1:numel(dirs)    
    % Find the title of this directory
    fidc = fopen('contents.txt','r+');
    prevdir = [];
    titletxt = fgetl(fidc);
    while ~strncmp(dirs{j},titletxt,numel(dirs{j}))
        prevdir = titletxt;
        titletxt = fgetl(fidc);
        if titletxt < 0, 
            error('CHEBFUN:examples:dirname', ...
                ['Unknown directory name "',dirs{j},'. Update contents.txt.']);
        end        
    end
    titletxt = titletxt(length(dirs{j})+2:end);
    
    % Add entry to examples/examples.html
    if ~strcmp(dirs{j},'temp')
        fprintf(fid0,['<li><a href="',dirs{j},'/',dirs{j},'.html">',titletxt,'</a>\n</li>\n\n']);
    end
    
    % Find the next and previous directories for breadcrumbs
    nextdir = fgetl(fidc);
    if isnumeric(nextdir)
        nextdir = [];
    else
        idx = strfind(nextdir,' ');
        nextdir = nextdir(1:idx-1);
    end
    if ~isempty(prevdir)
        idx = strfind(prevdir,' ');
        prevdir = prevdir(1:idx-1);
    end    
    fclose(fidc);
    if strcmp(nextdir,'temp'), nextdir = []; end
    if strcmp(dirs{j},'temp'), nextdir = []; prevdir = []; end
    
    % Move to the directory.
	cd(dirs{j})
    % Find all the m-files.
    dirlist = dir(fullfile(pwd,'*.m'));
    mfile = {dirlist.name};      
    
    % Make dirname/dirname.html
    fid = fopen([dirs{j},'.html'],'w+');
    % Write title
%     if shtml, fprintf(fid,['<div style="position:relative; left:-20px;">\n']);  end
    fprintf(fid,['                <h2>Chebfun Examples: ',titletxt,'</h2>\n']);
%     if shtml, fprintf(fid,'            </div>\n');  end
    if shtml
        % Make dirname/index.shtml
        make_shtml('index',dirs{j},[],titletxt,[],nextdir,prevdir);
    end
        
    % Get the ordering (we need to ignore A, AN, THE, ETC)
    desc = cell(numel(mfile),1); 
    for k = 1:numel(mfile)
        filename = mfile{k}(1:end-2);
               
        % Grab the file description.
        fidk = fopen([filename,'.m']);
        txt = fgetl(fidk);
        fclose(fidk);
        if txt < 1, continue, end % This mfile will be ignored
        if numel(txt) >1 && strcmp(txt(1:2),'%%')
            txt = txt(4:end);
        else
            txt = '     ';  % This mfile will be ignored
        end
        if strcmpi(txt(1:2),'A ')
            txt = txt(3:end);
        elseif strcmpi(txt(1:3),'AN ')
            txt = txt(4:end);            
        elseif strcmpi(txt(1:4),'THE ')
            txt = txt(5:end);
        end
        desc{k} = txt;
    end
    [desc indx] = sort(lower(desc));
    mfile = mfile(indx);
    
    % Loop over the files
    for k = 1:numel(mfile)
        cd(fullfile(examplesdir,dirs{j}));
        filename = mfile{k}(1:end-2);
               
        % Grab the file description (again).
        fidk = fopen([filename,'.m']);

        txt = fgetl(fidk); fclose(fidk);
        if txt < 1, continue, end % Ignore this file.
        if numel(txt) >1 && strcmp(txt(1:2),'%%')
            txt = txt(4:end);
        else
            continue % This mfile will be ignored
%             txt = [filename,'.m'];
        end
        if numel(txt) > 50
            idx = strfind(txt,':');
            if ~isempty(idx), txt = txt(1:idx(1)-1); end
        end
        fprintf(fid,['<span>',txt, '</span>     (']);
%         fprintf(fid,['<span style="text-transform:uppercase;">',txt, '</span>     (']);
        
        % Make dirname/html/filename.shtml
        if shtml
            if k < numel(mfile), next = mfile{k+1}(1:end-2); else next = []; end
            if k > 1, prev = mfile{k-1}(1:end-2); else prev = []; end
            make_shtml(filename,filename,'html',(txt),titletxt,next,prev);
        end
        
        if html || pdf
            % Make the example
            fprintf('Compiling %s/%s.m ...\n',dirs{j},filename);
            try 
                cd('../')
                make_examples(dirs{j},filename);
                cd(fullfile(examplesdir,dirs{j}));
            catch ME
                disp([dirs{j}, '/' ,filename ' CRASHED!'])
            end
        end

        % HTML link
        if shtml && exist(fullfile('html',[filename,'.shtml']),'file')
        % Link to dirname/html/filename.shtml
            fprintf(fid,'<a href="html/%s.shtml">html</a>, ',filename);
        elseif exist(fullfile('html',[filename,'.html']),'file')
        % Link to dirname/html/filename.html
            link = fullfile('html',[filename,'.html']);
            fprintf(fid,'<a href="%s">html</a>, ',link);
        end

        % Link to dirname/pdf/<filename>.pdf
        if shtml
            link = ['pdf/',filename,'.pdf'];
        else
            link = fullfile('pdf',[filename,'.pdf']);
        end
        fprintf(fid,'<a href="%s">PDF</a>, ',link);
        
        % Link to M-FILES %
        fprintf(fid,'<a href="%s.m">M-file</a>)\n',filename);
        fprintf(fid,'                <br/>\n\n');   

    end
    fclose(fid);

    fprintf([dirs{j},' published\n'])
    cd ..
end

% Open template
fid_et2 = fopen('templates/examples_template2.txt','r');
% Read data.
tmp = fread(fid_et2,inf,'*char');
fclose(fid_et2);
% Write
fprintf(fid0,' %s',tmp);
fclose(fid0);

% Upload to web server.
if shtml

    curdir = pwd;
    try
        cd(fullfile(webdir,'examples'))  
    catch
        return
    end
    
    for j = 1:numel(dirs)
        fprintf(['Uploading ',dirs{j},'. Please wait ... '])
        if ~exist(dirs{j},'dir'), mkdir(dirs{j}), end
        eval(['!cp ', fullfile(curdir,dirs{j},'*.m'),' ',dirs{j}])
        cd(dirs{j});
        eval(['!cp ', fullfile(curdir,dirs{j},[dirs{j},'.html']),' ',[dirs{j},'.html']])
        eval(['!cp ', fullfile(curdir,dirs{j},'index.shtml'),' ','index.shtml'])
        if exist(fullfile(curdir,dirs{j},'html'),'dir')
            eval(['!cp ', fullfile(curdir,dirs{j},'html','*.shtml'),' ','html'])
            if admin
                eval(['!chmod -Rf 775 ','html']); 
                eval(['!chgrp -Rf chebfun ','html/*']); 
                eval(['!chown -Rf hale ','html/*']); 
            end
        end
        cd ..
        if admin
            eval(['!chmod -f 775 *']);
            eval(['!chmod -Rf 775 ',dirs{j}]);
            eval(['!chgrp -Rf chebfun ',[dirs{j},'/*']]);
            eval(['!chown -Rf hale ',dirs{j}]);
        end  
        fprintf('Complete.\n')
    end
    cd(curdir)
    return
end

function make_shtml(file1,file2,dir,title,dirtitle,next,prev)
if nargin < 3, dir = []; end
% Create the files.
if ~isempty(dir)
    if ~exist(fullfile(dir),'dir')
        mkdir html
    end
    fid = fopen([dir,'/', file1,'.shtml'],'w+');
else
    fid = fopen([file1,'.shtml'],'w+');
end

% Open the templates.
fid1 = fopen('../templates/template1.txt','r');
% Template 3 has the footer contact info, whereas 2 doesn't.
if ~isempty(dir)
%     fid2 = fopen('../templates/template2.txt','r');
    fid2 = fopen('../templates/template3.txt','r');
else
    fid2 = fopen('../templates/template3.txt','r');
end

% Read their data.
tmp1 = fread(fid1,inf,'*char');
tmp2 = fread(fid2,inf,'*char');
fclose(fid1);
fclose(fid2);

% Write to file.
if ~isempty(dir)
    fprintf(fid,'<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">');
else
    fprintf(fid,'<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">');
    fprintf(fid,'\n<html xmlns="http://www.w3.org/1999/xhtml">\n');
end

% HEAD
fprintf(fid,'\n<head>\n');
% TITLE
fprintf(fid,'<title>%s</title>\n',(title)); % Lower case
% fprintf(fid,'<title>%s</title>\n',capitalize(title)); % Captialised

% META DATA
if ~isempty(dir)
    v = version; indx = strfind(v,'.'); v = v(1:indx(2)-1);
    fprintf(fid,'<meta name="generator" content="MATLAB %s">\n',v);
    fprintf(fid,'<link rel="schema.DC" href="http://purl.org/dc/elements/1.1/">\n');
    fprintf(fid,'<meta name="DC.date" content="%s">\n',datestr(now, 'yyyy-mm-dd'));
    fprintf(fid,'<meta name="DC.source" content="%s.m">\n',file1);
end

% TEMPLATE1
fwrite(fid, tmp1);

% BREADCRUMBS
if ~isempty(dir)
    % 2nd level
    fprintf(fid,'            <div id="breadcrumb">\n            <table id="bctable"><tbody><tr>\n');
    fprintf(fid,'            <td> > <a href="../../">examples</a> > <a href="../" class="lc">%s</a>',dirtitle);
    if ~isempty(prev)
        fprintf(fid, ' | <a href="%s.shtml">previous</a>',prev);
    else
        fprintf(fid, ' | previous');
    end
    if ~isempty(next)
        fprintf(fid, ' | <a href="%s.shtml">next</a>',next);
    else
        fprintf(fid, ' | next');
    end
    fprintf(fid,'</td>\n');
    fprintf(fid,'            <td style="text-align:right"> also available as: <a href="../pdf/%s.pdf">PDF</a> | <a href="../%s.m">M-file</a></td>\n',file1,file1);
    fprintf(fid,'            </tr></tbody></table>\n            </div>\n');
else
    % 1st level
    fprintf(fid,'            <div id="breadcrumb">\n            <table id="bctable"><tbody><tr>\n');
    fprintf(fid,'            <td> > <a href="../">examples</a>');
    if ~isempty(prev)
        fprintf(fid, ' | <a href="../%s/">previous</a>',prev);
    else
        fprintf(fid, ' | previous');
    end
    if ~isempty(next)
        fprintf(fid, ' | <a href="../%s/">next</a>',next);
    else
        fprintf(fid, ' | next');
    end
    fprintf(fid,'</td>\n            </tr></tbody></table>\n</div>\n');
end
% HTML INCLUDE
fprintf(fid,'            <!--#include virtual="%s.html" -->\n',file2);
% TEMPLATE 2
fwrite(fid, tmp2);
fclose(fid);


function mypublish(varargin)
close all
evalin('base','clear all');
chebfunpref('factory'), cheboppref('factory')
publish(varargin{:});
chebfunpref('factory'), cheboppref('factory')
close all

function txt = capitalize(txt)
txt = lower(txt);
txt(1) = upper(txt(1));
for k = 1:numel(txt)-1;
    tk = txt(k);
    if strcmp(tk,' ') || strcmp(tk,'-') || strcmp(tk,'(') || strcmp(tk,'[') || strcmp(tk,'/') || (~isempty(str2num(tk)) && isreal(str2num(tk)))
        txt(k+1) = upper(txt(k+1));
    end
end

txt = [txt '    '];

for k = 1:numel(txt)-4;
    tk = txt(k:k+4);
    if strncmp(tk,'A ',2)
        txt(k:k+1) = lower(tk(1:2));
    elseif strncmp(tk,'In ',3) || strncmp(tk,'Of ',3) || strncmp(tk,'At ',3) || strncmp(tk,'As ',3) || strncmp(tk,'An ',3)
        txt(k:k+2) = lower(tk(1:3));
    elseif strncmp(tk,'The ',4) || strncmp(tk,'And ',4) || strncmp(tk,'Its ',4) || strncmp(tk,'For ',4) || strncmp(tk,'Via ',4)
        txt(k:k+3) = lower(tk(1:4));
    elseif strcmp(tk,'With ')
        txt(k:k+4) = lower(tk);
    elseif strncmp(tk,'Ode',3) || strncmp(tk,'Pde',3) || strncmp(tk,'Bvp',3)
        txt(k:k+2) = upper(tk(1:3));
    end
end

txt = txt(1:end-4);
txt(1) = upper(txt(1));
% end

function make_tags

% Website directory
webdir = '/common/htdocs/www/maintainers/hale/chebfun/';
% Get tag info locally
examplesdir = pwd;

% The outputs
tagsfile = fullfile(examplesdir,'tags.php');
listfile = fullfile(examplesdir,'tags.list');

% Initialise the lists
tagslist = {}; filelist = {}; dirlist = {};
namelist = {}; tidynamelist = {};

% Find the relevent directories
alldirs = struct2cell(dir(fullfile(examplesdir)));
dirs = {}; k = 1;
for j = 3:size(alldirs,2)
    if alldirs{4,j} && ...
            ~strncmp(alldirs{1,j},'.',1) && ...
            ~strcmp(alldirs{1,j},'templates') && ...
            ~strcmp(alldirs{1,j},'old_examples') && ...
            ~strcmp(alldirs{1,j},'old') && ...
            ~strcmp(alldirs{1,j},'temp')
        dirs{k} = alldirs{1,j};
        k = k + 1;
    end
end

% Get all the files
for j = 1:size(dirs,2)
    tmp = struct2cell(dir(fullfile(examplesdir,dirs{j},'*.m')));
    tmp = tmp(1,:);
    filelist = [filelist tmp];
    tmp = repmat(dirs(j),1,numel(tmp));
    dirlist = [dirlist tmp];
end

for j = 1:numel(filelist);
    % Open the file
    fid = fopen(fullfile(examplesdir,dirlist{j},filelist{j}));
    
    % Get the titles
    txt = fgetl(fid);
    origtxt = txt;
    txt = upper(txt);
    if txt < 1, continue, end % This mfile will be ignored
    if numel(txt) >1 && strcmp(txt(1:2),'%%')
        txt = txt(4:end);
    else
        txt = '     ';  % This mfile will be ignored
    end
    if strcmpi(txt(1:2),'A ')
        txt = txt(3:end);
    elseif strcmpi(txt(1:3),'AN ')
        txt = txt(4:end);             
    elseif strcmpi(txt(1:4),'THE ')
        txt = txt(5:end);
    end
    tidynamelist{j} = txt;
    origtxt = origtxt(4:end);
    if numel(origtxt) > 50
        idx = strfind(origtxt,':');
        if ~isempty(idx), origtxt = origtxt(1:idx(1)-1); end
    end
    namelist{j} = origtxt;
        
    % Ignore these lines    
    for k = 2:5
        tmp = fgetl(fid);
    end
    tagline{j} = fgetl(fid);
    fclose(fid);
    
    % Seach for hashtags
    idx = strfind(tagline{j},'#');
    if isempty(idx)
        tagline{j} = '';
        continue
    end   
    idx2 = strfind(tagline{j},',');
    tags = {};
    for k = 1:numel(idx2)
        tags{k} = tagline{j}(idx(k)+1:idx2(k)-1);
    end
    tags{numel(idx)} = tagline{j}(idx(end)+1:end-1);
    tagline{j} = tagline{j}(idx(1):end-1);
    
    % Determine if new tag or existing
    for k = 1:numel(tags)
        if isempty(tagslist)
            idx = 0;
        else
            idx = strcmp(tags{k},tagslist(:,1));
        end
        if any(idx)
            idx = find(idx);
            tagslist{idx,2} = [tagslist{idx,2} j];
        else
            if isempty(tags{k})
                error('Error in tagging %s/%s',dirlist{j},filelist{j}), end
            tagslist{end+1,1} = tags{k};
            tagslist{end,2} = j;
        end
    end
end

% Sorting...
% Sort by number of times tagged
% [numtags idx] = sort(cellfun(@numel,tagslist(:,2)),'descend');
% Sort alphabetically
[ignored idx] = sort(tagslist(:,1));
% Do the sort
tagslist = tagslist(idx,:);

numtags = cellfun(@numel,tagslist(:,2));

[ignored idx2] = sort(tidynamelist);
% Write to the tags.list file
fid_list = fopen(listfile','w');
for k = 1:numel(tagline);
    j = idx2(k);
    % Get and clean the tags
    string = tagline{j};
    if isempty(string), continue, end
    string = [dirlist{j} '/' filelist{j}(1:end-2) ' ' strrep(string,',','') '\n'];
    fprintf(fid_list,string);
end
fclose(fid_list);


%% junk

% for j = 1:size(tagslist,1)
%     fprintf(fid,['<h2 id="' tagslist{j} '">#' tagslist{j} '</h2>\n']);
%     fprintf(fid,'<ul>\n');
%     kk = tagslist{j,2};
%     [ignored idx] = sort(tidynamelist(kk));
%     for k = kk(idx);
%         url = ['/chebfun/examples/' dirlist{k} '/html/' filelist{k}(1:end-2),'.shtml'];
%         fprintf(fid,['<li><a href="' url '">' namelist{k} '</a></li>\n']);
%     end
%     fprintf(fid,'</ul>\n\n');
% end

%%

% make some PHP/HTML magic
fid = fopen(tagsfile,'w+');

fidtmp = fopen([examplesdir '/templates/php_template1.txt'],'r');
while 1
    tline = fgetl(fidtmp);
    if ~ischar(tline), break, end
    fprintf(fid,tline);
    fprintf(fid,'\n');
end
fclose(fidtmp);

fprintf(fid,'<head>\n<title>Tag search</title>\n');
tmpdir = '/common/htdocs/www/maintainers/hale/chebfun/includes/';
fidtmp = fopen([tmpdir 'head.html'],'r');
while 1
    tline = fgetl(fidtmp);
    if ~ischar(tline), break, end
    if strfind(tline,'x-mathjax-config'), break, end
    fprintf(fid,tline);
    fprintf(fid,'\n');
end
fclose(fidtmp);
fprintf(fid,'</head>\n\n<body>\n<div id="headimg"></div>\n<div id="wrapper"> <!--WRAPPER-->');
fidtmp = fopen([tmpdir 'header.html'],'r');
while 1
    tline = fgetl(fidtmp);
    if ~ischar(tline), break, end
    fprintf(fid,tline);
    fprintf(fid,'\n');
end
fclose(fidtmp);
fidtmp = fopen([tmpdir 'mainmenu.html'],'r');
while 1
    tline = fgetl(fidtmp);
    if ~ischar(tline), break, end
    fprintf(fid,tline);
    fprintf(fid,'\n');
end
fclose(fidtmp);
fprintf(fid,'<h2 class="bigger">Tag search</h2>\n<div id="mycontent"> <!--CONTENT-->\n');

fprintf(fid,'<form action="tags.php" method="get">\n');
fprintf(fid,'search: <input type="text" name="query" />\n');
fprintf(fid,'<input type="submit" />\n</form>');

fidtmp = fopen([examplesdir '/templates/php_template2.txt'],'r');
while 1
    tline = fgetl(fidtmp);
    if ~ischar(tline), break, end
    fprintf(fid,tline);
    fprintf(fid,'\n');
end
fclose(fidtmp);


%%

fprintf(fid,['<h2 class="bigger" style="margin-left:-15pt; margin-top:10pt; margin-bottom:-10pt;">Tag cloud</h2>']);

% Create the cloud
maxfont = 500;
scl = ceil(maxfont/power(max(numtags),.666));
fsize = floor(max(scl*power(numtags,.66),50));
cloud = [];
for j = randperm(size(tagslist,1))
    if strcmp(tagslist{j},upper(tagslist{j})) && ...
       ~any(strcmpi(tagslist{j},{'FFT','2D','IVP'})) && ...
        isempty(str2num(tagslist{j}))
        continue
    end
    cloud = [cloud '<b style="font-size:' num2str(fsize(j)) '%%">',...
        '<a href="/chebfun/examples/tags.php?query=' tagslist{j} '">' tagslist{j} '</a></b>', ...
        '<b style="font-size:100%%"> </b>'];
%     '<a href="/chebfun/examples/tags.shtml#' tagslist{j} '">' tagslist{j} '</a></b>', ...
end
fprintf(fid,['<br/>\n\n<span style="width:400px">' cloud '</span>\n']);

fprintf(fid,['<h2 class="bigger" style="margin-left:-15pt; margin-top:10pt; margin-bottom:-10pt;">Function cloud</h2>']);
% Create the cloud
maxfont = 500;
scl = ceil(maxfont/power(max(numtags),.666));
fsize = floor(max(scl*power(numtags,.66),50));
cloud = [];
for j = randperm(size(tagslist,1))
    if ~strcmp(tagslist{j},upper(tagslist{j})) || ...
            any(strcmpi(tagslist{j},{'FFT','2D','IVP'})) || ...
            ~isempty(str2num(tagslist{j}))
            continue
    end
    cloud = [cloud '<b style="font-size:' num2str(fsize(j)) '%%">',...
        '<a href="/chebfun/examples/tags.php?query=' tagslist{j} '">' tagslist{j} '</a></b>', ...
        '<b style="font-size:100%%"> </b>'];
%     '<a href="/chebfun/examples/tags.shtml#' tagslist{j} '">' tagslist{j} '</a></b>', ...
end
fprintf(fid,['<br/>\n\n<span style="width:400px">' cloud '</span>\n']);

%%

fprintf(fid,'</div> <!--END CONTENT-->\n<div class="hr1"> <hr /> </div>\n');
fidtmp = fopen([tmpdir 'footer.html'],'r');
while 1
    tline = fgetl(fidtmp);
    if ~ischar(tline), break, end
    fprintf(fid,tline);
    fprintf(fid,'\n');
end
fclose(fidtmp);
fprintf(fid,'</div> <!--END WRAPPER-->\n<div id="footimg"></div>\n</body>\n</html>');
fclose(fid);

%%

% Upload
fprintf('Uploading tags...')
type(tagsfile)
type(listfile)
try
    eval(['!cp ' tagsfile ' ' webdir 'examples/'])
    eval(['!cp ' listfile ' ' webdir 'examples/'])
    currdir = pwd;
    eval(['!cd ' webdir 'examples/'])
    eval(['!chgrp -f chebfun ' webdir 'examples/tags.php'])
    eval(['!chgrp -f chebfun ' webdir 'examples/tags.list'])
    eval(['!cd ' currdir])
    fprintf('DONE!\n')
catch
    fprintf('FAILED!\n')
end

