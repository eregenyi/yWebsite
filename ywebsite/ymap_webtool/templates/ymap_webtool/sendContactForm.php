<!--send contact form-->
<!DOCTYPE html>
<html>
<body>
<?php

//VARIABLES
//$_POST is an array of variables passed to the current script via the HTTP POST method.
//Information sent from a form with the POST method is invisible to others (all names/values are embedded within the body of the HTTP request) and has no limits on the amount of information to send.
    $name = $_POST['name']; 
    $email = $_POST['email'];
    $message = $_POST['message'];
    $from = 'From: yMap'; 
    $to = 'paquayleila@hotmail.com'; 
    $subject = 'New message from a yMap user';
    $human = $_POST['human'];
			
    $body = "From: $name\n E-Mail: $email\n Message:\n $message";
				
    if ($_POST['submit'] && $human == '4') {				 
        if (mail ($to, $subject, $body, $from)) { 
	    echo '<p>Your message has been sent!</p>';
	} else { 
	    echo '<p>Something went wrong, go back and try again!</p>'; 
	} 
    } else if ($_POST['submit'] && $human != '4') {
	echo '<p>You answered the anti-spam question incorrectly!</p>';
    }
?>
</body>
</html>