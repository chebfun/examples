function new_examples(examples)
% Generates html for the 'new examples' on the Chebfun homepage.

% examples = {'geom','Ellipses'};

examplesdir = '/chebfun/examples/';
webdir = '/common/htdocs/www/maintainers/hale/chebfun/';

% color = {'#d9e3e5'
%          '#e4ecef'
%          '#eaefef'
%          '#eff2f2'};
     
color = {'#d9e3e5'
     '#e4ecef'
     '#eff2f2'};
     
if ischar(examples) || numel(examples) == 1 % Copy to the new exmaples page
    htmlfiles = dir(fullfile([webdir 'examples/new/'],'*.html'));
    newfilename = ['newexamples', int2str(numel(htmlfiles)+1), '.html'];
    cmd = ['!cp ',webdir, 'includes/newexamples.html ',webdir,'examples/new/',newfilename];
    eval(cmd)
    fid = fopen('newexamples.html','w+');
    
    text = fileread([webdir,'examples/new/',newfilename]);
    text = strrep(text,'<h3>New Examples:</h3>','');
    text = strrep(text,'exampleImage',['exampleImage' int2str(numel(htmlfiles)+1)]);
    text = strrep(text,'imageLink',['imageLink' int2str(numel(htmlfiles)+1)]);
    fprintf(fid,text);
    fprintf(fid,'\n\n</br>\n\n');
    j = 1;
    
    htmlnames = {htmlfiles.name};
    htmlnames = strrep(htmlnames,'newexamples','');
    htmlnames = strrep(htmlnames,'.html','');
    for k = 1:numel(htmlnames)
        htmlnames{k} = str2num(htmlnames{k});
        %if isempty(htmlnames{k}), htmlnames{k} = NaN; end
    end
    [ignored idx] = sort([htmlnames{:}]);
    htmlfiles = htmlfiles([1 idx+1]);
    
    for k = numel(htmlfiles):-1:1
        if strcmp(htmlfiles(k).name,'newexamples.html'), continue, end
        name = strrep(htmlfiles(k).name,'newexamples',''); name = strrep(name,'.html','');

        text = fileread([webdir,'examples/new/',htmlfiles(k).name]);
        text = strrep(text,'<h3>New Examples:</h3>','');
        text = strrep(text,'exampleImage',['exampleImage' name]);
        text = strrep(text,'imageLink',['imageLink' name]);
        
        j = j+1;
        if ~mod(j,2) % flip left and right
            text = strrep(text,'style="float:right;','tempOX26Le');
            text = strrep(text,'style="float:left;','style="float:right;');
            text = strrep(text,'tempOX26Le','style="float:left;');
            text = strrep(text,'thumbnail','thumbnailleft');
            text = strrep(text,'<span><br/>','<span style="margin-left:-300px;">br/>');
        end
        
        fprintf(fid,text);
        fprintf(fid,'\n\n</br>\n\n');
        
    end
    fclose(fid);
    cmd = ['!mv newexamples.html ',webdir,'examples/new/'];
    eval(cmd)
    return
end

fid = fopen('newexamples.html','w+');
dir1 = examples{1}; file1 = examples{2};

newline = '\n'     ;
% str1 = '        <h3>New Examples:</h3>';
str1 = '';
str2 = '        <div style="background:#fff; color:#000; height:160px;">';
str3 = '            <div style="float:left; background:#fff; width:605px; padding:1px 0px; height:160px">';
str = [str1,newline,str2,newline,str3];
fprintf(fid,str);

defaultimage = [];
k = 0;
while ~isempty(examples)
    k = k+1;
    dirk = examples{1};
    file = examples{2};
    if numel(examples) > 2 && ~ischar(examples{3})
        imageno = examples{3};
        if isempty(imageno), imageno = 1; end
        imageno = int2str(imageno); 
        if length(imageno) == 1, imageno = ['0' imageno]; end
        examples(1:3) = [];
    else
        imageno = '01';
        examples(1:2) = [];
    end
    if isempty(defaultimage), defaultimage = imageno; end

    filedir = [examplesdir, dirk, '/html/'];
    STYLE = ['style="background:',color{k},';"'];
    tmp = [dirk,'/',file,'.m'];
    fidk = fopen([dirk,'/',file,'.m']);
    try
        title = fgetl(fidk);
    catch ME
        [dirk,'/',file,'.m']
        rethrow(ME)
    end
    title = title(4:end);
    nameanddate = fgetl(fidk);
    nameanddate = nameanddate(3:end);
    fclose(fidk);
    
    ws = '                ';
    str1 = '<b class="rtop2" style="margin-top:3px"><b class="r1" STYLE></b><b class="r2" STYLE></b><b class="r3" STYLE></b><b class="r4" STYLE></b></b>';
    str1 = strrep(str1,'STYLE',STYLE);
    fprintf(fid,[newline, ws, str1]);
    
    str2 = '<div style="margin:0px; padding:4px; padding-left:20px; background:';
    fprintf(fid,[newline,ws,str2,color{k},';">']);
    
    str3 = ['<a href="',filedir,file,'.shtml" ',...
            'onmouseover="document.images[''exampleImage''].src=''',filedir,file,'_',imageno,'.png'';',...
                         'document.images[''exampleImageB''].src=''',filedir,file,'_',imageno,'.png'';',...
                          'document.links[''imageLink''].href=''',filedir,file,'.shtml''"',...
                         '>',title,'</a>'];
    str3b = ['<p style="font-size:11px; margin:0px; padding:0px;">',nameanddate,'</p>'];
    fprintf(fid,[newline, ws,'   ',str3,newline, ws,'   ',str3b,newline,ws,'</div>']);
    
    str4 = '<b class="rbottom2"><b class="r4" STYLE></b><b class="r3" STYLE></b><b class="r2" STYLE></b><b class="r1" STYLE></b></b>';
    str4 = strrep(str4,'STYLE',STYLE);
    fprintf(fid,[newline, ws, str4]);
end

str = {};
str{1} = '            </div>';
str{2} = '            <div style="float:right; width:220px; padding-left:0px; padding-right:0px; margin-top:-10px;">';
str{3} = '                <b class="rtop2" style="margin-top:5px"><b class="r1"></b><b class="r2"></b><b class="r3"></b><b class="r4"></b></b>';
str{4} = ['                <div style="float:right; background:',color{1},'; width:190; height:160px; padding-left:10px; padding-right:10px; align:center;">'];
str{5} = ['                    <a class="thumbnail" name="imageLink" href="',examplesdir,dir1,'/html/',file1,'.shtml">',...
                              '<img src="',examplesdir,dir1,'/html/',file1,'_',defaultimage,'.png" height="150px" width="200px" name="exampleImage" style="margin-top:5px;" border="0px">',...
                              '<span style="margin-left:-300px;"><br/><img src="',examplesdir,dir1,'/html/',file1,'_',defaultimage,'.png" name="exampleImageB"/></span></a>'];
str{6} = '                </div>';
str{7} = '                <b class="rbottom2" style="margin-top:160px;"><b class="r4"></b><b class="r3"></b><b class="r2"></b><b class="r1"></b></b>';
str{8} = '            </div>';
str{9} = '      </div>';

for k = 1:numel(str)
    fprintf(fid,[newline,str{k}]);
end

fclose(fid);

cmd = ['!mv newexamples.html ',webdir,'includes/'];
eval(cmd)

