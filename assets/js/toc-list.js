$(document).ready(function() {
    var tocList = $('#toc-list');
    $('.post-content').find('h2, h3').each(function() {
      var tag = $(this).prop('tagName').toLowerCase();
      var text = $(this).text();
      var id = text.toLowerCase().replace(/ /g, '-').replace(/[^\w-]+/g, '');
      
      // 给每个标题添加id
      $(this).attr('id', id);

      // 根据标题层级生成大纲
      var li = $('<li></li>').append('<a href="#' + id + '">' + text + '</a>');
      if(tag === 'h3') {
        li.css('margin-left', '20px'); // 子标题缩进
      }
      tocList.append(li);
    });
  });