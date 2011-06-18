#!/usr/bin/env python
'''Downloads the Madagascar wiki as static html
Wiki pages are downloaded every time. Main_Page is moved to index.html
Html files in local dir are deleted if not present on the wiki, in order
to avoid accumulation of dead wiki pages. To save bandwidth cost and
download time, images and linked documents are downloaded only if a
local copy is not present.'''

# Copyright (C) 2008, 2009 Ioan Vlad
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

from urllib import urlopen
import os, copy, string, sys

try: # Give precedence to local version
    import ivlad
except: # Use distributed version
    import rsf.user.ivlad as ivlad

###############################################################################

def before(str,sep):
    '''
    return part before the separator
    '''
    return str.split(sep)[0]

def after(str,sep):
    '''
    return part after the separator
    '''
    strs = str.split(sep)+['']
    return strs[1]

###############################################################################

def generate_hard_coded_paths(outdir):
    h = {}

    # The wiki should be found at wikirooturl, downloaded to wikiroot_local
    h['hosting_domain'] = 'http://ahay.org'
    h['wiki_basenm']    = 'wiki'
    h['wiki_url']       = h['hosting_domain'] + '/' + h['wiki_basenm']+'/'
    h['wiki_local']     = outdir

    # Images should be in some subdirs of imgdir_url, downloaded to imgdir_local
    h['imgdir_basenm'] = 'images'
    h['imgdir_url']    = h['wiki_url'] + h['imgdir_basenm']
    h['imgdir_local']  = os.path.join(h['wiki_local'], h['imgdir_basenm'])

    # Docs should be at docdir_url (and subdirs), downloaded to docdir_local
    h['docdir_basenm'] = 'docs'
    h['nonwiki_root']  = h['hosting_domain'] + '/wikilocal'
    h['docdir_url']    = h['nonwiki_root'] + '/' + h['docdir_basenm']
    h['docdir_local']  = os.path.join(h['wiki_local'], h['docdir_basenm'])

    # Stylesheet to be found at wikiroot/css
    h['css_basenm'] = 'commonPrint.css'
    h['css_relpath']= 'skins/common/' + h['css_basenm']
    h['path_from_wikiroot'] = h['wiki_basenm'] + '/' + h['css_relpath']
    h['css_url']    = h['wiki_url'] + h['css_relpath']
    h['css_local']  = os.path.join(h['wiki_local'], h['css_basenm'])

    # PHP argument in the wiki URL
    h['printable'] = '&printable=yes'

    # Extension for output files
    h['ext'] = '.html'

    return h

###############################################################################

def read_page_section(hcp, pagenm, just_before_section, just_after_section):
    'Reads a page from wiki and extracts a section from it'
    f = urlopen(hcp['wiki_url'] + pagenm + hcp['printable'])
    s = f.read()
    f.close()
    s2 = after(s,just_before_section)
    return before(s2,just_after_section)

###############################################################################

def get_page_list(hcp, verb):
    'Returns a list with pages to download from the wiki'

    ivlad.msg('Downloading and parsing page list...', verb)

    just_before_pagelist = '<hr />'
    just_before_pagelist += \
    '<table style="background: inherit;" border="0" width="100%">'
    just_before_pagelist += '<tr><td width="33%">'

    just_after_pagelist='</table><div class="printfooter">'

    s = read_page_section(hcp, 
                          'Special:AllPages', 
                          just_before_pagelist, 
                          just_after_pagelist)
    list1 = s.split('<td width="33%">')
    p0 = '<a href="'
    p1 = '" title="'
    p2 = '/%s/' % hcp['wiki_basenm']
    list2 = map(lambda x: after(before(x.lstrip(p0),p1),p2), list1)
    ivlad.msg('...done', verb)
    return sorted(list2)

###############################################################################

def clean_and_mkdir_local(hcp, pagelist):
    'Create local dir, or clean up pages removed from wiki'
    d = hcp['wiki_local']
    if os.path.isdir(d):
        files = sorted(os.listdir(d))
        for f in files:
            (f_nm, f_ext) = os.path.splitext(f)
            if f_ext == hcp['ext']:
                if f_nm not in pagelist:
	            os.remove(os.path.join(d,f))
    else: # local dir does not exist
        os.mkdir(d)

    if not os.path.isdir(hcp['imgdir_local']):
        os.mkdir(hcp['imgdir_local'])
    if not os.path.isdir(hcp['docdir_local']):
        os.mkdir(hcp['docdir_local'])

###############################################################################

def get_wiki_image_dict(hcp, verb):
    'Returns a dictionary with image_name:path/inside/img/url/dir'

    ivlad.msg('Downloading and parsing image list...', verb)

    pg = 'Special:ImageList?limit=500&ilsearch=&title=Special%3AImageList'
    just_before_imglist = '<th>Description</th>\n</tr></thead><tbody>'
    just_after_imglist = '''</tr>
</tbody></table>
<br />
<table class="imagelist_nav TablePager_nav" align="center" cellpadding="3">'''
    just_after_imglist += '<tr>'
    img_table = read_page_section(hcp, 
                                  pg, 
                                  just_before_imglist, 
                                  just_after_imglist)
    list1 = img_table.split('<tr>\n<td class="TablePager_col_img_timestamp">')
    imgdict = {}
    p1 = '<td class="TablePager_col_img_name"><a href="/%s/Image:' % \
        hcp['wiki_basenm']
    p2 = '">file</a>)</td>'
    p3 = '(<a href="/%s/%s/' % (hcp['wiki_basenm'],hcp['imgdir_basenm'])

    for item in list1:
	    s = after(before(after(item,p1),p2),p3)
	    if s != '':
	        (filepath, filename) = os.path.split(s)
	        imgdict[filename] = filepath

    ivlad.msg('...done', verb)
    return imgdict

###############################################################################

def inspect_imgs(hcp, wiki_img_dict):
    'Returns dict with urls to download, targets'
    wid = copy.copy(wiki_img_dict)
    local_img_list   = sorted(os.listdir(hcp['imgdir_local']))
    wiki_img_list    = sorted(wiki_img_dict.keys())
    for img in local_img_list:
        if img in wiki_img_list:
	    # Do not download image again. Imgs practically never change, and
	    # bandwidth costs Madagascar community members after-tax money
            del wid[img]
	else: # image has been removed from the wiki. remove locally too
	    os.path.remove(os.path.join(hcp['imgdir_local'],img))
    imgs_to_download = {}
    for i in wid:
        url    = hcp['imgdir_url'] + '/' + wid[i] + '/' + i
        target = os.path.join(hcp['imgdir_local'], i)
        imgs_to_download[url] = target
    return imgs_to_download

###############################################################################

def build_wiki_img_repl_dict(hcp, wiki_img_dict):
    'Dictionary with image hyperlinks replacements in wiki pages'
    wiki_img_repl_dict = {}
    widk = wiki_img_dict.keys()
    s1 = 'src="/%s/%s/' % (hcp['wiki_basenm'],hcp['imgdir_basenm'])
    s2 = 'src="./%s/' % hcp['imgdir_basenm']
    for img in widk:
        to_be_replaced = s1 + wiki_img_dict[img] + '/' + img
	replacement    = s2 + img
        wiki_img_repl_dict[to_be_replaced] = replacement
    return wiki_img_repl_dict

###############################################################################

def download_file(url, filename, verb, whatisfile, textmode=False):
    'Downloads a file from the web'
    format = 'w'
    if not textmode:
        format += 'b'
    ivlad.msg('Downloading %s: %s' % (whatisfile, url), verb)
    urlhandle = urlopen(url)
    contents  = urlhandle.read()
    outhandle = open(filename, format)
    outhandle.write(contents)
    outhandle.close()

###############################################################################

def filter_page(s, hcp, pagelist, wiki_img_repl_dict):
    's is a string containing a html page'
    # links to docs
    s=s.replace(hcp['nonwiki_root'],'.')
    # link to style sheet
    css_str_to_replace  = 'href="/%s?164"' % hcp['path_from_wikiroot']
    css_str_replacement = 'href="%s"'      % hcp['css_basenm']
    s = s.replace(css_str_to_replace,css_str_replacement)
    # links to other pages
    for pagenm in pagelist:
        str_to_replace  = '<a href="/' + hcp['wiki_basenm'] + '/' + pagenm
        replacement_str = '<a href="' + pagenm + hcp['ext']
        s = s.replace(str_to_replace,replacement_str)
    # links to images
    link_beg = '<a href="/%s/Image:' % hcp['wiki_basenm']
    for imgpath in wiki_img_repl_dict.keys():
        s = s.replace(imgpath,wiki_img_repl_dict[imgpath])
        img = os.path.split(imgpath)[1]
        imgtitle = img.replace('_',' ')
        lbeg = link_beg + img + '" class="image" title="'
        s = s.replace(lbeg+         imgtitle        +'">', '')
        s = s.replace(lbeg+'Image:'+imgtitle        +'">', '')
        s = s.replace(lbeg+'Image:'+imgtitle.lower()+'">', '')
        s = s.replace(lbeg+'Image:'+img             +'">', '')
        s = s.replace(lbeg+'Image:'+img.lower()     +'">', '')
    return s

###############################################################################

def download_pages(hcp, pagelist, wiki_img_repl_dict):
    for page in pagelist:
        ivlad.msg('Downloading and parsing: ' + page)
        page_url    = hcp['wiki_url'] + page + hcp['printable']
        if page == 'Main_Page':
            local_basename = 'index'
        else:
            local_basename = page
        if '/' in local_basename:
            # Temporary fix -- this will break some links...
            local_basename = string.replace(local_basename, '/', '')
        page_html = os.path.join(hcp['wiki_local'], local_basename+hcp['ext'])
        page_handle = urlopen(page_url)
        s = page_handle.read()
        page_handle.close()
        s = filter_page(s,hcp,pagelist,wiki_img_repl_dict)
        f = open(page_html, 'w')
        f.write(s)
        f.close()

###############################################################################

def download_imgs(imgs_to_download, wikiroot, imgdir):
    for imgfile in sorted(imgs_to_download.keys()):
        ivlad.msg('Downloading: ' + imgfile)
        img_url = wikiroot + img_dir_basename + '/' + \
                  imgs_to_download[imgfile] + '/' + imgfile
        img_file = os.path.join(imgdir,imgfile)
        ivlad.exe('curl ' + img_url + '>' + img_file)

###############################################################################

def inspect_doc_dir(docdir_url, docdir_local, docs_to_download):
    just_before_section = '/">Parent Directory</a>'
    just_after_section = '<hr></pre>'
    f = urlopen(docdir_url)
    s1 = f.read()
    f.close()
    s2 = after(s1,just_before_section)
    s3 = before(s2,just_after_section)
    list1 = s3.split('<a href="')
    list2 = map(lambda x: before(x,'">').strip(),list1)
    list2.remove('-')
    for entry in list2:
        local_node = os.path.join(docdir_local, entry)
        www_node = docdir_url
        if docdir_url[-1] != '/':
             www_node += '/'
        www_node += entry
        if entry[-1] == '/': # it's a directory
            if not os.path.isdir(local_node):
                os.mkdir(local_node)
            inspect_doc_dir(www_node, local_node, docs_to_download)
        else:
            if not os.path.isfile(local_node):
                docs_to_download[www_node] = local_node
###############################################################################

def inspect_docs(hcp, verb):
    'generates with docs to download (url:target)'
    ivlad.msg('Downloading and parsing document list...', verb)
    docs_to_download = {}
    inspect_doc_dir(hcp['docdir_url'], hcp['docdir_local'], docs_to_download)
    ivlad.msg('...done', verb)
    return docs_to_download

###############################################################################

def main(par):

    verb = par.bool('verb', False) # verbosity flag
    outdir = par.string('outdir','.')

    hcp = generate_hard_coded_paths(outdir)

    pagelist = get_page_list(hcp, verb)
    clean_and_mkdir_local(hcp, pagelist)

    wiki_img_dict      = get_wiki_image_dict(     hcp, verb)
    imgs_to_download   = inspect_imgs(            hcp, wiki_img_dict)
    wiki_img_repl_dict = build_wiki_img_repl_dict(hcp, wiki_img_dict)

    docs_to_download = inspect_docs(hcp, verb)

    download_pages(hcp, pagelist, wiki_img_repl_dict)

    for i in imgs_to_download:
        download_file(i, imgs_to_download[i], verb, 'image')

    for d in docs_to_download:
        download_file(d, docs_to_download[d], verb, 'document')

    if not os.path.isfile(hcp['css_local']):
        download_file(hcp['css_url'], 
                      hcp['css_local'], 
                      verb, 
                      'stylesheet', 
                      textmode=True)

    return ivlad.unix_success

###############################################################################

if __name__ == '__main__':
    ivlad.run(main)
