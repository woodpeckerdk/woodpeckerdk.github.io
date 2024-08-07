$(document).ready(function() {
    var toc = $('#toc'); // 获取目录容器
    var headings = $('.content h2, .content h3').filter(function() {
    return $(this).text().trim() !== ''; // 过滤掉空标题
    });
 
    if (headings.length > 0) {
    var tocList = $('<ul></ul>'); // 创建目录列表
    var currentH2 = null;
 
    headings.each(function(index, heading) {
        var headingText = $(heading).text();
        var headingId = 'heading-' + index;
        $(heading).attr('id', headingId); // 为每个标题添加唯一ID
 
        var listItem = $('<li></li>').append(
        $('<a></a>').attr('href', '#' + headingId).text(headingText)
        );
 
        if ($(heading).is('h2')) {
            listItem.addClass('toc-h2');
            tocList.append(listItem); // 将h2目录项添加到列表中
            currentH2 = listItem;
        } else if ($(heading).is('h3') && currentH2) {
            listItem.addClass('toc-h3'); // 为h3目录项添加样式类
            if (!currentH2.find('ul').length) {
                currentH2.append($('<ul></ul>'));
            }
            currentH2.find('ul').append(listItem); // 将h3目录项添加到对应的h2子列表中
        }
    });
 
    toc.append(tocList); // 将目录列表添加到目录容器中
    }
});
