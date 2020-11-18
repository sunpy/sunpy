#!/bin/bash
timestamp=$(date '+%F %H:%M').strip()
message=f"Azure Pipelines Build for {$BRANCH} {$STATUS} at {timestamp} (https://dev.azure.com/sunpy/sunpy/_build/results?buildId={$BUILDID})"
colour="#00ff00" if $STATUS == "Succeeded" else "#ff0000"
formatted=f"Azure Pipelines Build for {$BRANCH} <font color=\"{colour}\">{$STATUS}</font> at <a href='https://dev.azure.com/sunpy/sunpy/_build/results?buildId={$BUILDID}'>{timestamp}</a>"
url=f"{$HOMESERVER}/_matrix/client/r0/rooms/{$ROOMID}/send/m.room.message"

echo @(message)
echo @(formatted)

http --ignore-stdin POST @(url) f"Authorization: Bearer {$ACCESS_TOKEN}" msgtype=m.notice f"body={message}" format=org.matrix.custom.html f"formatted_body={formatted}"
