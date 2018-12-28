#-*- mode: ruby ; coding: utf-8 -*-
#= Generate HTML from RD
#
# Authors:: Youhei SASAKI
# License:: Expat (a.k.a. MIT/X11)
#+DATE: 2013-09-03 20:02:31

require 'rubygems'
require 'rake'
require 'yaml'
require 'rd/rdfmt'
require 'rd/rd2html-ext-lib'
require 'RMagick'
require 'image_size'
require 'liquid'
require 'pp'
require 'pathname'
require 'nkf'


#== set BASE Directory
#
BASE = File.expand_path(File.dirname(__FILE__))
#
#== load config.yaml
#
config = YAML::load_file("#{BASE}/config.yml")

if config["base_path"] == nil
  config["base_path"] = "/library/dcmodel/htmltools/Webgen"
end
liquid = Hash.new
liquid["site"] = config

######################################################################
# source files
if config["recursive"] == true
  ALL_SRC_FILES = FileList["#{BASE}/**/*.rd"]
elsif config["src"] != (nil|false)
  ALL_SRC_FILES = FileList["#{BASE}/*.rd"]
  unless config["src"]
    config["src"].each do |s|
      ALL_SRC_FILES.include("#{BASE}/#{s}/*.rd")
    end
  end
else
  ALL_SRC_FILES = FileList["#{BASE}/*.rd"]
end
######################################################################
# all target files
GENERATED_FILES = []; ALL_SRC_FILES.each do |fpath|
  basename = fpath.gsub("\.rd", '')
  GENERATED_FILES.push basename + ".htm.ja" if config["lang"]["ja"]
  GENERATED_FILES.push basename + ".htm.en" if config["lang"]["en"]
end

######################################################################
desc "Create/Update html files"
task :default => GENERATED_FILES

######################################################################
# FileList["#{BASE}/**/*_thumbnail*"] if config["thumbnail"]

desc "=> clobber"
task :clean => :clobber
desc "=> clobber"
task :distclean => :clobber
desc "clean all generated files"
task :clobber do
  GENERATED_FILES.each {|fn| rm_r fn rescue nil}
end

#= 拡張子ルール
#
rule '.ja' => "%{\.htm$,.rd}X"  do |t|
  puts "Update #{t} from #{t.source}"
  RD2Webgen.rdswap "#{t}", "#{t.source}", lang="ja", liquid
end
rule '.en' => "%{\.htm$,.rd}X"  do |t|
  puts "Update #{t} from #{t.source}"
  RD2Webgen.rdswap "#{t}", "#{t.source}", lang="en", liquid
end

######################################################################
private

module RD2Webgen

  #== HEADER, FOOTER
  #
  # 外部から読み込めるようにも設定するが, 決め打ち様に
  #
  HEADER = <<-HEAD_EOF
<!DOCTYPE html>
<html lang="{{ site.lang }}">
  <head>
    <meta charset="UTF-8" />
    <meta name="author" content="{{ site.author }}">
    <title>{{ site.title }}: {{ site.page_title }}</title>
    <meta name="viewport" content="width=device-width initial-scale=1.0">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <link rel="stylesheet" href="{{ site.base_path }}/assets/bootstrap/css/bootstrap.min.css" type="text/css" medial="all">
    <link rel="stylesheet" href="{{ site.base_path }}/assets/bootstrap/css/bootstrap-theme.min.css" type="text/css" medial="all">
    <!-- HTML5 Shim and Respond.js IE8 support of HTML5 elements and media queries -->
    <!-- WARNING: Respond.js doesn't work if you view the page via file:// -->
    <!--[if lt IE 9]>
       <script src="{{ site.base_path }}/assets/html5shiv/html5shiv.js"></script>
       <script src="{{ site.base_path }}/assets/respond.js/respond.min.js"></script>
    <![endif]-->
    {% if site.css %}
    <link href="{{ site.css }}" rel="stylesheet">
    {% endif %}
    <script src="{{ site.base_path }}/assets/jquery/jquery.min.js"></script>
    <script src="{{ site.base_path }}/assets/bootstrap/js/bootstrap.min.js"></script>
    {% if site.additional_head %}
    {{ site.additional_head }}
    {% endif %}
  <style type="text/css">
   .navbar-nav > li {
     font-size: 10pt;
   }
   .navbar-nav > li > a {
     padding-left: 2px;
     padding-right: 2px;
     padding-top: 16px;
     padding-bottom: 12px;
   }
   .navbar-nav > li > a:before {
     content: " / ";
     color: #ccc;
   }
   .navbar-brand > a {
   }
   .navbar-form {
     margin-right: 5px;
     margin-left: 5px;
   }
  </style>
  </head>
  <body>
  <nav class="navbar navbar-default navbar-static-top" role="navigation">
  <div class="navbar-header">
    <button type="button" class="navbar-toggle" data-toggle="collapse" data-target="#header01">
      <span class="sr-only">Toggle navigation</span>
      <span class="icon-bar"></span>
      <span class="icon-bar"></span>
      <span class="icon-bar"></span>
    </button>
    <a class="navbar-brand" href="{{ site.url }}">{{ site.title }}</a>
  </div>
  <div class="collapse navbar-collapse" id="header01">
    <ul class="nav navbar-nav">
      {% if site.lang == 'ja' %}
      {% for link in site.link_ja %}
      <li><a href="{{ link.url }}">{{ link.title }}</a></li>
      {% endfor %}
      <li><a href="SIGEN.htm">SIGEN</a></li>
      {% else %}
      {% for link in site.link_en %}
      <li><a href="{{ link.url }}">{{ link.title }}</a></li>
      {% endfor %}
      {% endif %}
    </ul>
    <ul class="nav navbar-nav navbar-right">
    <li>
    {% if site.lang == 'ja' %}
    <a href="{{ site.page_url || replace:'.htm.ja','.htm.en' }}">English</a>
    {% else %}
    <a href="{{ site.page_url || replace:'.htm.en','.htm.ja' }}">日本語</a>
    {% endif %}
    </li>
    <li>
    <form class="navbar-form" role="search"
                 action="http://google.com/search"
                 method="get"
                 onsubmit="$('#siteseach').val($location).attr('host')+$(location).attr('patname')">
        <input type="hidden" name="q"  value="site:{{ site.url | remove_first:'http://' }}">
        <input type="hidden" name="ie" value="UTF-8">
        <input type="hidden" name="oe" value="UTF-8">
        <input type="hidden" name="sitesearch" id="sitesearch">
        {% if site.lang == 'ja' %}
        <input type="text" name="q"  size="6" class="form-control" placeholder="検索">
        {% else %}
        <input type="text" name="q"  size="6" class="form-control" placeholder="search">
        {% endif %}
    </form>
    </li>
    </ul>
    </div><!-- /.navbar-collapse -->
    </nav>
    <div class="container main">
HEAD_EOF
  FOOTER = <<-FOOT_EOF
  </div> <!-- /.container main -->
  <div class="container main">
  <footer>
    <p>
    &copy; {{ site.author }}: <img src="{{ site.email }}" height=16pt alt="contact address" style="display-style: inline; vertical-align: middle;" /> <br />
    Lastupdate: {{ site.mtime }}, Since {{ site.since }}
    </p>
  </footer>
  </div> <!-- /.container main -->
  </body>
</html>

FOOT_EOF
  #== GC
  #
  # RMagick(ImageMagick) は結構メモリを喰うので, 明示的に GC しないと危険
  #
  def self.run_gc
    fDisabled = GC.enable
    GC.start
    GC.disable if fDisabled
  end

  #== make thumbnail
  #
  # RD ファイル内の IMG リンクから画像のファイル名を取り出して,
  # thumbnail を作成する. 今の所 PNG 決め打ちであることに注意されたい.
  #
  def self.make_rmagick_thumbnail(file, input)
    @file = file
    @input = input
    scale = 0.3
    _src_base_path = File.dirname(@file)
    @input.split("\n").each do |l|
      if l =~ /\(\(<(.*)\|"IMG:(.*)">\)\)$/
        _alt_name = $1; _img_url = $2.gsub("./","")
        unless _img_url =~/_thumbnail/
          _orig_img_path = File.join(_src_base_path, _img_url)
          _thumb_img_name = "/" + File.basename(_img_url).gsub(/\.(.+?)$/, '_thumbnail.\1')
          _thumb_img_path = File.join(_src_base_path, File.dirname(_img_url), _thumb_img_name)
          unless File.exists?(_thumb_img_path)
            puts "Create thumbnail #{_thumb_img_path} from #{_orig_img_path}"
            img = Magick::Image.read(_orig_img_path).first
            img.thumbnail!(scale)
            img.format = "PNG"
            File.open(_thumb_img_path, "w+"){|f|
              f.puts img.to_blob
            }
            run_gc()
          end
        end
      end
    end
  end

  #
  #== Rakefile の置き場所からの相対パスを計算する
  #
  def self.calc_relative_path(path, from)
    @path = Pathname.new(path)
    @from = Pathname.new(from)
    return @path.relative_path_from(@from).to_s
  end

  #
  #== Liquid テンプレートのコンパイル
  #
  def self.compile_liquid(liquid)
    @liquid = liquid
    header = Liquid::Template.parse(HEADER).render(@liquid)
    footer = Liquid::Template.parse(FOOTER).render(@liquid)
    return header, footer
  end

  #== RD ファイルの処理
  #
  # 日英併記された RD ファイルから HTML ファイルを生成する
  #
  def self.rdswap(dest, file, lang, liquid, debug=nil)
    @file = file
    @dest = dest
    @contents = open(@file)
    @lang = lang
    @liquid = liquid
    @source = NKF.nkf('-w', @contents.read).gsub(/\r\n/, "\n")
    @contents.close
    # 必要に応じてサムネイルを作成
    # make_rmagick_thumbnail(@file, @source) if @liquid["site"]["thumbnail"]
    # RD ファイルの処理
    tree = RD::RDTree.new(@source, [RD::RDTree.tmp_dir], nil)
    visitor = RD::RD2HtmlDennouClubVisitor.new
    # visitor.install_ref_extension
    # visitor.install_native_inline
    # tree.filter["RT"] = RD::RT_FILTER
    # visitor.include_suffix.push("RT")
    tree.filter["HTML"] = RD::INCLUDE_FILTER
    visitor.include_suffix.push("HTML")
    if @lang == "ja"
      tree.filter["JA"] = RD::RD_FILTER
      visitor.include_suffix.push("JA")
      tree.filter["HTMLJA"] = RD::INCLUDE_FILTER
      visitor.include_suffix.push("HTMLJA")
      @liquid["site"]["lang"] = "ja"
    elsif @lang == "en"
      tree.filter["EN"] = RD::RD_FILTER
      visitor.include_suffix.push("EN")
      tree.filter["HTMLEN"] = RD::INCLUDE_FILTER
      visitor.include_suffix.push("HTMLEN")
      @liquid["site"]["lang"] = "en"
    end
    begin
      result = visitor.visit(tree.parse)
    rescue
      puts "result: #{@dest} is empty."
      @liquid["site"]["page_title"] = "404 Not Found"
      _page_url = File.basename(@dest)
      @liquid["site"]["page_url"] = _page_url
      @liquid["site"]["mtime"] = File.mtime(@file).strftime("%Y/%m/%d")
      header, footer = compile_liquid(@liquid)
      open(@dest, "w+"){|f|
        f.puts header
        f.puts "<p>now under construction</p>"
        f.puts footer
      }
    else
      result.each_line do |l|
        if l =~/<h1>(.*)<\/h1>/
          _page_title = $1.gsub(/<a\sname=\".*\">(.*)<\/a>/, '\1')
          @liquid["site"]["page_title"] = _page_title
        end
      end
      _page_url = File.basename(@dest)
      @liquid["site"]["page_url"] = _page_url
      src_base_path = File.dirname(@file)
      rel_path = calc_relative_path(BASE, File.dirname(@dest))
      # if @liquid["site"]["publish_base"]
      #   @liquid["site"]["assets"] = @liquid["site"]["publish_base"] + "assets"
      # else
      #   @liquid["site"]["assets"] = rel_path + "/assets"
      # end
      @liquid["site"]["mtime"] = File.mtime(@file).strftime("%Y/%m/%d")
      # Liquid テンプレートのコンパイル
      header, footer = compile_liquid(@liquid)
      # 出力
      open(@dest, "w+"){|f|
        f.puts header
        f.puts result
        f.puts footer
      }
    end
    # 一次ファイルの削除
    Dir.glob("#{RD::RDTree.tmp_dir}/rdtmp.#{$$}.*.*").each do |i|
      File.delete(i)
    end
  end


end

######################################################################
# RD extension
module RD

  class RD2HtmlDennouClubVisitor < RD::RD2HTMLExtVisitor

    # module RD::RD2HTMLExtVisitor::RefExtension
      # alias :old_ref_ext_IMG :ref_ext_IMG

      # private
      # def ref_ext_IMG(element, label, content)
      #   return nil unless /^IMG:(.+)$/i =~ label
      #   file = $1
      #   if @liquid["site"]["thumbnail"]
      #     thumbnail = file.gsub(/(.*)\.(.+?)$/, '\1_thumbnail.\2')
      #     ret = %Q[<a href="#{file}">]
      #     ret << %Q[<img src="#{thumbnail}" alt="#{content}" image_size_should_be_replaced />]
      #     ret << %Q[</a>]
      #   else
      #     ret = %Q[<img src="#{file}" alt="#{content}" class="img-responsive"/>]
      #   end
      # end
    # end

    def visit(tree)
      install_ref_extension
      install_native_inline
      super
    end

    # def install_native_inline
    #   extend NativeInline
    # end

    # remove "id"
    def apply_to_Headline(element, title)
      anchor = get_anchor(element)
      label = hyphen_escape(element.label)
      title = title.join("")
      %Q[<h#{element.level}><a name="#{anchor}">#{title}</a></h#{element.level}>]
    end # apply_to_Headline

    def apply_to_DescListItem(element, term, description)
      anchor = get_anchor(element.term)
      label = hyphen_escape(element.label)
      term = term.join("")
      if description.empty?
        %Q[<dt><a name="#{label}">#{term}</a></dt>]
      else
        %Q[<dt><a name="#{label}">#{term}</a></dt>\n] +
          %Q[<dd>\n#{description.join("\n").chomp}\n</dd>]
      end
    end # apply_to_DescListItem

    def apply_to_DocumentElement(element, contents)
      content = contents.join("\n")
      foottext = make_foottext
      snippet = "<section class=\"contents\">\n#{content}\n</section>\n"
      if foottext
        foottext = foottext.gsub(/\A<hr \/>/, '')
        foottext = "<h2>Footnotes</h2>\n#{foottext}"
        snippet << "<section class=\"footnotes\">\n#{foottext}\n</section>\n"
      end
      snippet = "<div id=\"main\">\n#{snippet}</div>"
      snippet
    end # apply_to_DocumentElement

    def apply_to_DocumentElement(element, contents)
      return content = contents.join("\n")
    end

  end

  # RT_FILTER = RD::Filter.new(:target) do |inn, out|
  #   out.print("=\begin RT\n")
  #   inn.each do |line|
  #     out.print(line)
  #   end
  #   out.print("\n=\end RT\n")
  # end #RT_FILTER

  HTML_FILTER = RD::Filter.new(:target) do |inn,out|
    out.print("=begin HTML\n")
    inn.each do |line|
      out.print(line)
    end
    out.print("=\n\end HTML\n")
  end #HTML_FILTER

end
