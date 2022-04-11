#!/bin/bash
timestamp=$(date '+%F %H:%M').strip()
message=f"GitHub Actions Build for {$BRANCH} {$STATUS} at {timestamp} (https://github.com/sunpy/sunpy/actions/runs/{$RUNID})"
colour="#00ff00" if $STATUS == "Succeeded" else "#ff0000"
formatted=f"GitHub Actions Build for {$BRANCH} <font color=\"{colour}\">{$STATUS}</font> at <a href='https://github.com/sunpy/sunpy/actions/runs/{$RUNID}'>{timestamp}</a>"
url=f"{$HOMESERVER}/_matrix/client/r0/rooms/{$ROOMID}/send/m.room.message"

echo @(message)
echo @(formatted)

http --ignore-stdin POST @(url) f"Authorization: Bearer {$ACCESS_TOKEN}" msgtype=m.notice f"body={message}" format=org.matrix.custom.html f"formatted_body={formatted}"
