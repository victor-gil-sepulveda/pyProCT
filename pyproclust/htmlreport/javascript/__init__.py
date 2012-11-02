
# From http://jasalguero.com/ledld/development/web/expandable-list/ !!!

page_header_chunk = """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
        <title>Demo Expandable list</title>
        <link rel="stylesheet" href="css/style.css" type="text/css" media="screen, projection">
    <link rel="stylesheet" href="css/jquery.rotateTableCellContent.css"/>
        
    <script type="text/javascript" src="js/jquery-1.4.2.min.js"></script>
        <script type="text/javascript" src="js/scripts.js"></script>
        <script type="text/javascript" src="js/jquery.rotateTableCellContent.js"></script>
    
    <script type="text/javascript">
        $(document).ready(function(){
            $('#summary_table').rotateTableCellContent();
        });

    </script>

    </head>
    <body>
"""

page_footer_chunk = """
    </body>
</html>
"""